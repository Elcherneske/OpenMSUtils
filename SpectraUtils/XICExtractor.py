from math import ceil
from typing import List
import pickle
from dataclasses import dataclass
from tqdm import tqdm
import pandas as pd
import multiprocessing
from multiprocessing import shared_memory
import numpy as np
from .MSUtils import SpectraObject, MSReader

# 配置日志
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# to accelerate the parallel processing, we use a global variable to store the peaks
shared_peaks_buf_name = None
shared_peaks_offsets = None

def preprocess_fn(preprocessor):
    return preprocessor.preprocess()

@dataclass
class XICResult:
    """XIC 结果类"""
    rt_array: list[float]  # 保留时间数组
    intensity_array: list[float]  # 强度数组
    mz: float  # 目标质荷比
    ppm_array: list[float]  # PPM 误差

class XICExtractor:
    """XIC 提取器类"""
    def __init__(
        self, 
        ppm_tolerance: float = 25.0, 
        rt_bin_size: float = 1.0, 
        num_threads: int = 1, 
        min_scans: int = 5,
        peak_boundary: float = 0.2,
        mode: str = 'rt_range',
    ):
        """
        初始化 XIC 提取器
        
        参数:
            mzml_file: mzML 文件路径
            ppm_tolerance: 质量容差 (ppm)
            rt_bin_size: 保留时间分箱大小
            num_threads: 并行处理的线程数s
        """
        self.ppm_tolerance = ppm_tolerance
        self.rt_bin_size = rt_bin_size
        self.num_threads = num_threads
        self.mode = mode
        self.min_scans = min_scans
        self.peak_boundary = peak_boundary

        # 是否调整XIC的起止点，默认不调整
        if self.mode == 'scan_window':
            self.adjust_xic = True
        else:
            self.adjust_xic = False
        
    def load_mzml(self, mzml_file: str):
        """加载 mzML 文件并转换为 MSObject 列表"""

        def format_ms_clusters(spectra_objects: list[SpectraObject], rt_bin_size: float):
            """
            将 MSObject 列表转换为簇格式, 并记录每个cluster的rt
            
            参数:
                ms_objects: MSObject 列表, 按rt排序
            """

            global shared_peaks_buf_name, shared_peaks_offsets

            # 按rt排序
            spectra_objects.sort(key=lambda x: x.retention_time)

            # create shared memory for peaks
            total_peaks_num = sum(len(spectra_obj.peaks) for spectra_obj in spectra_objects)
            peaks_array = np.empty(shape=(total_peaks_num, 2), dtype=np.float32)
            current_offset = 0
            offsets = []
            for spectra_obj in spectra_objects:
                n = len(spectra_obj.peaks)
                if n == 0:
                    offsets.append((current_offset * 8, current_offset * 8))
                    continue
                peaks_array[current_offset:current_offset + n, :] = spectra_obj.peaks
                start_byte = current_offset * 8
                end_byte = (current_offset + n) * 8
                offsets.append((start_byte, end_byte))
                current_offset += n
            
            shared_peaks_buf = shared_memory.SharedMemory(create=True, size=peaks_array.nbytes)
            assert peaks_array.nbytes == current_offset * 8
            shared_peaks_buf.buf[:peaks_array.nbytes] = peaks_array.tobytes()
            shared_peaks_buf_name = shared_peaks_buf.name
            shared_peaks_offsets = offsets

            rt_indices = [(-1, -1)] * (int(spectra_objects[-1].retention_time / rt_bin_size) + 1)
            ms_clusters = []

            current_cluster = None
            for index, spectra_obj in enumerate(spectra_objects):
                if spectra_obj.level == 1:
                    if not current_cluster:
                        current_cluster = {'rt': spectra_obj.retention_time, 'ms1': index, 'ms2': []}
                    else:
                        ms_clusters.append(current_cluster)
                        rt_index = int(current_cluster['rt'] / rt_bin_size)
                        if rt_indices[rt_index][0] == -1 or rt_indices[rt_index][1] == -1:
                            rt_indices[rt_index] = (len(ms_clusters) - 1, len(ms_clusters) - 1)
                        rt_indices[rt_index] = (rt_indices[rt_index][0], len(ms_clusters) - 1)
                        current_cluster = {'rt': spectra_obj.retention_time, 'ms1': index, 'ms2': []}
                else:
                    current_cluster['ms2'].append({'mz_min': spectra_obj.precursor_window[0], 'mz_max': spectra_obj.precursor_window[1], 'index': index, 'rt':spectra_obj.retention_time})
            
            if current_cluster:
                ms_clusters.append(current_cluster)
                rt_index = int(current_cluster['rt'] / rt_bin_size)
                if rt_indices[rt_index][0] == -1 or rt_indices[rt_index][1] == -1:
                    rt_indices[rt_index] = (len(ms_clusters) - 1, len(ms_clusters) - 1)
                rt_indices[rt_index] = (rt_indices[rt_index][0], len(ms_clusters) - 1)
            
            return {
                'ms_clusters': ms_clusters,
                'rt_indices': rt_indices,
            }

        reader = MSReader(thread_num=self.num_threads)
        spectra_objects = reader.read_to_spectra_objects(mzml_file)
        
        if not spectra_objects:
            raise ValueError("未能从 mzML 文件中读取到谱图数据")

        ms_clusters_infos = format_ms_clusters(spectra_objects, self.rt_bin_size)
        return ms_clusters_infos
    
    def extract_xics(self, mzml_file: str, df: pd.DataFrame) -> List[tuple[List[XICResult], List[XICResult]]]:
        """提取XIC"""
        logging.info(f"load mzml file {mzml_file} with {self.num_threads} threads...")
        ms_clusters_infos = self.load_mzml(mzml_file)
        rt_indices = ms_clusters_infos['rt_indices']
        ms_clusters = ms_clusters_infos['ms_clusters']

        # 预处理XIC条目
        logging.info(f"preprocess xic entries with {self.num_threads} threads...")
        xic_entries = []
        if self.num_threads > 1:
            chunk_size = max(1, ceil(len(df) / self.num_threads))
            chunks = [df.iloc[i:i + chunk_size] for i in range(0, len(df), chunk_size)]

            # 创建预处理器
            preprocessors = []
            for i, chunk in enumerate(chunks):
                preprocessor = RangePreprocessor(
                    df=chunk,
                    rt_indices=rt_indices,
                    ms_clusters=ms_clusters,
                    rt_bin_size=self.rt_bin_size,
                    mode=self.mode
                )
                preprocessors.append(preprocessor)
            
            with multiprocessing.Pool(processes=self.num_threads) as pool:
                results = list(tqdm(
                    pool.imap(preprocess_fn, preprocessors),
                    total=len(preprocessors),
                    desc=f"预处理XIC条目"
                ))
                # 收集处理结果
                for result in results:
                    if result:
                        xic_entries.extend(result)
        else:
            preprocessor = RangePreprocessor(
                df=df,
                rt_indices=rt_indices,
                ms_clusters=ms_clusters,
                rt_bin_size=self.rt_bin_size,
                mode=self.mode
            )
            xic_entries = preprocessor.preprocess()

        # 提取XIC
        logging.info(f"extract xics with {self.num_threads} threads...")
        xics = []
        if self.num_threads > 1:            
            # 将数据划分为多个块
            chunk_size = max(1, ceil(len(xic_entries) / self.num_threads))
            chunks = [xic_entries[i:i + chunk_size] for i in range(0, len(xic_entries), chunk_size)]

            preprocessors = []
            for i, chunk in enumerate(chunks):
                preprocessor = ExtractPreprocessor(
                    xic_entries=chunk,
                    ppm_tolerance=self.ppm_tolerance,
                    min_scans=self.min_scans,
                    peak_boundary=self.peak_boundary,
                    adjust_xic=self.adjust_xic
                )
                preprocessors.append(preprocessor)
            
            # 创建进程池并并行处理
            with multiprocessing.Pool(processes=self.num_threads) as pool:
                results = list(tqdm(
                    pool.imap(preprocess_fn, preprocessors),
                    total=len(chunks),
                    desc=f"提取XIC (使用{self.num_threads}个进程)"
                ))
                
                # 合并结果
                for result in results:
                    xics.extend(result)
        else:
            # 单线程处理
            preprocessor = ExtractPreprocessor(
                xic_entries=xic_entries,
                ppm_tolerance=self.ppm_tolerance,
                min_scans=self.min_scans,
                peak_boundary=self.peak_boundary,
                adjust_xic=self.adjust_xic
            )
            xics = preprocessor.preprocess()
        
        global shared_peaks_buf_name, shared_peaks_offsets
        # 释放两个shared memory
        if shared_peaks_buf_name is not None:
            try:
                from multiprocessing import shared_memory
                shm = shared_memory.SharedMemory(name=shared_peaks_buf_name)
                shm.close()
                shm.unlink()
                shared_peaks_buf_name = None
            except Exception as e:
                logging.warning(f"Error releasing shared_peaks_buf_name '{shared_peaks_buf_name}': {e}")
        if shared_peaks_offsets is not None:
            shared_peaks_offsets = None

        return xics

class RangePreprocessor:
    def __init__(
            self, 
            df: pd.DataFrame,
            rt_indices: list[tuple[int, int]], 
            ms_clusters: list[dict], 
            rt_bin_size: float,
            mode: str
        ):
        self.df = df
        self.rt_indices = rt_indices
        self.ms_clusters = ms_clusters
        self.rt_bin_size = rt_bin_size
        self.mode = mode

    def preprocess(self):
        def calculate_isotope_mzs(mz: float, charge: int, num_isotopes: int = 4) -> List[float]:
            """计算同位素峰的质荷比"""
            isotope_mzs = [mz]
            for i in range(1, num_isotopes):
                # 碳同位素间隔约为 1.0033 Da
                isotope_mz = mz + (i * 1.00335483507) / charge
                isotope_mzs.append(isotope_mz)
            return isotope_mzs
        
        xic_entries = []
        for index, row in self.df.iterrows():
            precursor_mzs = calculate_isotope_mzs(float(row["precursor_mz"]), int(row["charge"]))
            fragment_mzs = [float(mz) for mz in row["fragment_mz"].split(",")]

            if self.mode == 'rt_range':
                assert 'rt_start' in row and 'rt_stop' in row, "rt_range mode requires rt_start and rt_stop columns"
                rt_start = float(row["rt_start"]) * 60
                rt_stop = float(row["rt_stop"]) * 60
                ms_cluster = self._get_ms_scan_by_rt_range((rt_start, rt_stop), precursor_mzs[0])
            elif self.mode == 'scan_window' or self.mode == 'fix_scan_window':
                assert 'rt' in row and 'scan_window_num' in row, "scan_window mode requires rt and scan_window_num columns"
                rt = float(row["rt"]) * 60
                scan_window_num = int(row["scan_window_num"])
                ms_cluster = self._get_ms_scan_by_scan_window(rt, scan_window_num, precursor_mzs[0])
            else:
                raise ValueError(f"Invalid mode: {self.mode}")
            
            xic_entries.append(
                {
                    "precursor_mzs": precursor_mzs,
                    "fragment_mzs": fragment_mzs,
                    "ms_cluster": ms_cluster
                }
            )
        return xic_entries

    def _get_ms_scan_by_scan_window(self, rt: float, scan_window_num: int, precursor_mz: float) -> list[dict]:
        """
        获取指定保留时间范围内的峰
        
        参数:
            rt: 保留时间
            scan_window_num: 扫描窗口大小, scan_window_num = upper - lower - 1 = (2 * radius)
        """

        # get search index of ms_clusters 
        rt_bin = int(rt / self.rt_bin_size)
        rt_low_index = -1
        rt_high_index = -1

        bin_index = rt_bin
        while self.rt_indices[bin_index][0] == -1 or self.rt_indices[bin_index][1] == -1:
            bin_index -= 1
            if bin_index < 0:
                break
        if bin_index < 0:
            rt_low_index = 0
        else:
            rt_low_index = self.rt_indices[bin_index][0]

        bin_index = rt_bin
        while self.rt_indices[bin_index][0] == -1 or self.rt_indices[bin_index][1] == -1:
            bin_index += 1
            if bin_index > len(self.rt_indices) - 1:
                break
        if bin_index > len(self.rt_indices) - 1:
            rt_high_index = len(self.rt_indices) - 1
        else:
            rt_high_index = self.rt_indices[bin_index][1]
        
        # get the closest ms_cluster
        selected_index = rt_low_index
        min_rt_diff = float('inf')
        for index in range(rt_low_index, rt_high_index + 1):
            rt_diff = abs(self.ms_clusters[index]['rt'] - rt)
            if rt_diff < min_rt_diff:
                min_rt_diff = rt_diff
                selected_index = index
        
        low_index = selected_index - scan_window_num//2
        high_index = selected_index + scan_window_num//2
        if low_index < 0:
            low_index = 0
        if high_index >= len(self.ms_clusters):
            high_index = len(self.ms_clusters) - 1
        
        clusters = self.ms_clusters[low_index:high_index + 1]
        greedy_index = 0
        results = []
        for cluster in clusters:
            result = {'rt': cluster['rt'], 'ms1': cluster['ms1'], 'ms2': None}
            filtered_ms2 = None
            if 0 <= greedy_index < len(cluster['ms2']) and cluster['ms2'][greedy_index]['mz_min'] <= precursor_mz <= cluster['ms2'][greedy_index]['mz_max']:
                filtered_ms2 = cluster['ms2'][greedy_index]['index']
            else:
                for index, ms2 in enumerate(cluster['ms2']):
                    if ms2['mz_min'] <= precursor_mz <= ms2['mz_max']:
                        filtered_ms2 = ms2['index']
                        greedy_index = index
                        break
            result['ms2'] = filtered_ms2
            results.append(result)
        return results
    
    def _get_ms_scan_by_rt_range(self, rt_range: tuple[float, float], precursor_mz: float) -> list[dict]:
        """
        获取指定保留时间范围内的峰
        
        参数:
            rt_range: (rt_min, rt_max) 元组
            
        返回:
            峰列表，每个元素为 (mz, intensity) 元组
        """
        rt_low, rt_high = rt_range
        bin_low = int(rt_low / self.rt_bin_size)
        bin_high = int(rt_high / self.rt_bin_size)
        low_index = -1
        high_index = -1
        for bin_index in range(bin_low, bin_high + 1):
            if self.rt_indices[bin_index][0] != -1 and self.rt_indices[bin_index][1] != -1:
                if low_index == -1:
                    low_index = self.rt_indices[bin_index][0]
                    break
        
        for bin_index in range(bin_high, bin_low - 1, -1):
            if self.rt_indices[bin_index][0] != -1 and self.rt_indices[bin_index][1] != -1:
                if high_index == -1:
                    high_index = self.rt_indices[bin_index][1]
                    break
        
        if low_index == -1 or high_index == -1:
            clusters = self.ms_clusters
        else:
            if low_index > 0:
                low_index -= 1 # 有可能 low_index-1 的 ms1 rt不在范围，但是 ms2 rt在范围
            clusters = self.ms_clusters[low_index:high_index + 1]
        
        # filter clusters by rt range
        # clusters = [cluster for cluster in clusters if cluster['rt'] >= rt_low and cluster['rt'] <= rt_high]
        greedy_index = 0
        results = []
        for cluster in clusters:
            result = {'ms1': None, 'ms2': None}

            # select ms1
            selected_ms1 = None
            if rt_low <= cluster['rt'] <= rt_high:
                selected_ms1 = {'index': cluster['ms1'], 'rt': cluster['rt']}
            result['ms1'] = selected_ms1

            # select ms2
            selected_ms2 = None
            if (0 <= greedy_index < len(cluster['ms2']) 
                and cluster['ms2'][greedy_index]['mz_min'] <= precursor_mz <= cluster['ms2'][greedy_index]['mz_max']
                and rt_low <= cluster['ms2'][greedy_index]['rt'] <= rt_high
            ):
                selected_ms2 = {'index': cluster['ms2'][greedy_index]['index'], 'rt': cluster['ms2'][greedy_index]['rt']}
            else:
                for index, ms2 in enumerate(cluster['ms2']):
                    if (ms2['mz_min'] <= precursor_mz <= ms2['mz_max']
                        and rt_low <= ms2['rt'] <= rt_high
                    ):
                        selected_ms2 = {'index': ms2['index'], 'rt': ms2['rt']}
                        greedy_index = index
                        break
            result['ms2'] = selected_ms2
            if selected_ms1 is not None or selected_ms2 is not None:
                results.append(result)
            else:
                continue
        return results
    
class ExtractPreprocessor:
    def __init__(
            self, 
            xic_entries: list[dict],
            ppm_tolerance: float,
            min_scans: int,
            peak_boundary: float,
            adjust_xic: bool
        ):
        self.xic_entries = xic_entries
        self.ppm_tolerance = ppm_tolerance
        self.min_scans = min_scans
        self.peak_boundary = peak_boundary
        self.adjust_xic = adjust_xic
    
    def preprocess(self) -> list[tuple[List[XICResult], List[XICResult]]]: 
        """从mz簇中提取XIC"""
        def extract_xic_from_peaks(peaks_list: list[np.ndarray], mz: float, ppm_tolerance: float) -> tuple[List[float], List[float]]:
            """从峰列表中提取XIC"""
            def binary_search(peaks: np.ndarray, mz: float, ppm_tolerance: float) -> tuple[int, float]:
                """使用二分法找到最后一个比mz小的peak对应的index"""
                left, right = 0, len(peaks) - 1
                target_idx = 0  # 默认为列表长度，表示所有元素都小于mz
                
                while left <= right:
                    mid = (left + right) // 2
                    if peaks[mid][0] <= mz:
                        target_idx = mid
                        left = mid + 1
                    else:
                        right = mid - 1
                
                # 从target_idx开始，找最接近的peak
                min_ppm = float('inf')
                index = target_idx
                for i in range(target_idx, len(peaks)):
                    ppm = abs(peaks[i][0] - mz) / mz * 1e6
                    if ppm < min_ppm:
                        min_ppm = ppm
                        index = i
                    if (peaks[i][0] - mz)/mz * 1e6 > ppm_tolerance:
                        break
                
                if min_ppm < ppm_tolerance:
                    return peaks[index][1], min_ppm
                else:
                    return 0, ppm_tolerance

            intensity_list = []
            ppm_list = []
            
            for peaks in peaks_list:
                if len(peaks) == 0:
                    intensity_list.append(0)
                    ppm_list.append(ppm_tolerance)
                else:
                    intensity, ppm = binary_search(peaks, mz, ppm_tolerance)
                    intensity_list.append(intensity)
                    ppm_list.append(ppm)
            
            return intensity_list, ppm_list

        global shared_peaks_buf_name, shared_peaks_offsets

        peaks_buf = shared_memory.SharedMemory(name=shared_peaks_buf_name)

        xics = []
        for xic_entry in self.xic_entries:
            precursor_mzs = xic_entry["precursor_mzs"]
            fragment_mzs = xic_entry["fragment_mzs"]
            ms_cluster = xic_entry["ms_cluster"]

            #precursor_mz的XIC
            precursor_xics = []
            selected_clusters = [cluster for cluster in ms_cluster if cluster['ms1'] is not None]
            ms1_rt_list = [cluster['ms1']['rt'] for cluster in selected_clusters]
            ms1_specs_offsets = [shared_peaks_offsets[cluster['ms1']['index']] for cluster in selected_clusters]
            ms1_specs = [np.ndarray(shape=((offset[1] - offset[0])//8, 2), buffer=peaks_buf.buf[offset[0]:offset[1]], dtype=np.float32) for offset in ms1_specs_offsets]
            for precursor_mz in precursor_mzs:
                precursor_intensity_list, precursor_ppm_list = extract_xic_from_peaks(ms1_specs, precursor_mz, self.ppm_tolerance)
                precursor_xics.append(XICResult(ms1_rt_list, precursor_intensity_list, precursor_mz, precursor_ppm_list))
            
            #fragment_mz的XIC
            fragment_xics = []
            selected_clusters = [cluster for cluster in ms_cluster if cluster['ms2'] is not None]
            ms2_rt_list = [cluster['ms2']['rt'] for cluster in selected_clusters]
            ms2_specs_offsets = [shared_peaks_offsets[cluster['ms2']['index']] for cluster in selected_clusters]
            ms2_specs = [np.ndarray(shape=((offset[1] - offset[0])//8, 2), buffer=peaks_buf.buf[offset[0]:offset[1]], dtype=np.float32) for offset in ms2_specs_offsets]
            for fragment_mz in fragment_mzs:
                if all(len(ms2_spec) == 0 for ms2_spec in ms2_specs):
                    fragment_xics.append(XICResult(ms2_rt_list, [0] * len(ms2_rt_list), fragment_mz, [self.ppm_tolerance] * len(ms2_rt_list)))
                else:
                    fragment_intensity_list, fragment_ppm_list = extract_xic_from_peaks(ms2_specs, fragment_mz, self.ppm_tolerance)
                    fragment_xics.append(XICResult(ms2_rt_list, fragment_intensity_list, fragment_mz, fragment_ppm_list))

            if self.adjust_xic:
                precursor_xics, fragment_xics = self._adjust_xic(precursor_xics, fragment_xics)

            xics.append((precursor_xics, fragment_xics))
        return xics

    def _adjust_xic(self, precursor_xics: List[XICResult], fragment_xics: List[XICResult]) -> tuple[List[XICResult], List[XICResult]]:
        """调整XIC的起止点
        
        参数:
            precursor_xics: 前体离子XIC列表
            fragment_xics: 碎片离子XIC列表
            
        返回:
            调整后的前体离子和碎片离子XIC列表
        """
        def find_peak_boundaries(xic: XICResult) -> tuple[int, int]:
            """找到单个XIC的峰边界
            
            参数:
                xic: XIC结果对象
                
            返回:
                (start_idx, stop_idx): 峰起止点的索引
            """
            # 对intensity进行平滑处理
            intensity_array = xic.intensity_array
            smoothed_intensity = [0.0] * len(intensity_array)
            
            # 平滑处理
            smoothed_intensity[0] = (2.0 / 3.0) * intensity_array[0] + (1.0 / 3.0) * intensity_array[1]
            smoothed_intensity[-1] = (2.0 / 3.0) * intensity_array[-1] + (1.0 / 3.0) * intensity_array[-2]
            for i in range(1, len(intensity_array) - 1):
                smoothed_intensity[i] = 0.5 * intensity_array[i] + 0.25 * (intensity_array[i - 1] + intensity_array[i + 1])
            
            # 使用中间点作为最大强度点
            mid_idx = len(smoothed_intensity) // 2
            max_intensity = smoothed_intensity[mid_idx]
            
            # 设置边界阈值
            boundary = max_intensity * self.peak_boundary
            
            # 向左搜索起点
            start_idx = mid_idx
            valley = max_intensity
            valley_pos = mid_idx
            
            for i in range(mid_idx - 1, -1, -1):
                if smoothed_intensity[i] < valley:
                    valley = smoothed_intensity[i]
                    valley_pos = i
                elif valley < max_intensity / 3.0 and valley < smoothed_intensity[i] / 2.0:
                    start_idx = valley_pos
                    break
                if smoothed_intensity[i] < boundary:
                    break
            
            # 向右搜索终点
            stop_idx = mid_idx
            valley = max_intensity
            valley_pos = mid_idx
            
            for i in range(mid_idx + 1, len(smoothed_intensity)):
                if smoothed_intensity[i] < valley:
                    valley = smoothed_intensity[i]
                    valley_pos = i
                elif valley < max_intensity / 3.0 and valley < smoothed_intensity[i] / 2.0:
                    stop_idx = valley_pos
                    break
                if smoothed_intensity[i] < boundary:
                    break
            
            return start_idx, stop_idx
        
        # 收集所有XIC的边界
        all_boundaries = []
        for xic in precursor_xics + fragment_xics:
            if len(xic.rt_array) >= self.min_scans:
                start_idx, stop_idx = find_peak_boundaries(xic)
                if stop_idx - len(xic.rt_array)//2 < self.min_scans//2 or len(xic.rt_array)//2 - start_idx < self.min_scans//2:
                    continue
                all_boundaries.append((start_idx, stop_idx))
        
        # 如果有足够的XIC，计算最大边界
        if len(all_boundaries) > 0:
            max_start = max([b[0] for b in all_boundaries])
            max_stop = min([b[1] for b in all_boundaries])
        else:
            max_start = len(xic.rt_array)//2 - self.min_scans//2
            max_stop = len(xic.rt_array)//2 + self.min_scans//2
        
        # 调整所有XIC
        adjusted_precursor_xics = []
        for xic in precursor_xics:
            if len(xic.rt_array) > 0:
                start_idx = max(0, max_start)
                stop_idx = min(len(xic.rt_array) - 1, max_stop)
                adjusted_xic = XICResult(
                    rt_array=xic.rt_array[start_idx:stop_idx + 1],
                    intensity_array=xic.intensity_array[start_idx:stop_idx + 1],
                    mz=xic.mz,
                    ppm_array=xic.ppm_array[start_idx:stop_idx + 1]
                )
                adjusted_precursor_xics.append(adjusted_xic)
            else:
                adjusted_precursor_xics.append(xic)
        
        adjusted_fragment_xics = []
        for xic in fragment_xics:
            if len(xic.rt_array) > 0:
                start_idx = max(0, max_start)
                stop_idx = min(len(xic.rt_array) - 1, max_stop)
                adjusted_xic = XICResult(
                    rt_array=xic.rt_array[start_idx:stop_idx + 1],
                    intensity_array=xic.intensity_array[start_idx:stop_idx + 1],
                    mz=xic.mz,
                    ppm_array=xic.ppm_array[start_idx:stop_idx + 1]
                )
                adjusted_fragment_xics.append(adjusted_xic)
            else:
                adjusted_fragment_xics.append(xic)
        
        return adjusted_precursor_xics, adjusted_fragment_xics