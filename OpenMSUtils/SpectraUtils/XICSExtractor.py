from math import ceil
import numpy as np
from typing import List
from dataclasses import dataclass
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from .MZMLUtils import MZMLReader
from .MSObject import MSObject

@dataclass
class XICResult:
    """XIC 结果类"""
    rt_array: list[float]  # 保留时间数组
    intensity_array: list[float]  # 强度数组
    mz: float  # 目标质荷比
    ppm_array: list[float]  # PPM 误差

class XICSExtractor:
    """XIC 提取器类"""
    def __init__(
        self, 
        mzml_file: str, 
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
            num_threads: 并行处理的线程数
        """
        self.mzml_file = mzml_file
        self.ppm_tolerance = ppm_tolerance
        self.rt_bin_size = rt_bin_size
        self.num_threads = num_threads
        self.mode = mode
        self.min_scans = min_scans
        self.peak_boundary = peak_boundary

        # 根据mode设置是否调整XIC
        if self.mode == 'rt_range':
            self.adjust_xic = False
        elif self.mode == 'scan_window':
            self.adjust_xic = True
        else:
            self.adjust_xic = False
        
        self.ms_clusters = []
        self.rt_indices = {}
        
    def load_mzml(self):
        """加载 mzML 文件并转换为 MSObject 列表"""
        print(f"正在读取 mzML 文件: {self.mzml_file}")
        reader = MZMLReader()
        ms_objects = reader.read_to_msobjects(self.mzml_file, parallel=True)
        
        if not ms_objects:
            raise ValueError("未能从 mzML 文件中读取到谱图数据")
        
        print(f"共读取 {len(ms_objects)} 个谱图, 构建ms1,ms2谱图列表")  
        self._format_ms_clusters(ms_objects)
    
    def extract_xics(self, df: pd.DataFrame) -> List[tuple[List[XICResult], List[XICResult]]]:
        """提取XIC"""
        xic_entries = self._format_xic_entries(df)
        # 提取XIC
        xics = []
        if self.num_threads > 1:
            import multiprocessing as mp
            
            # 将数据划分为多个块
            chunk_size = max(1, ceil(len(xic_entries) / self.num_threads))
            chunks = [xic_entries[i:i + chunk_size] for i in range(0, len(xic_entries), chunk_size)]
            
            # 创建进程池并并行处理
            with mp.Pool(processes=self.num_threads) as pool:
                results = list(tqdm(
                    pool.imap(self._batch_extract_xic, chunks),
                    total=len(chunks),
                    desc=f"提取XIC (使用{self.num_threads}个进程)"
                ))
                
                # 合并结果
                for result in results:
                    xics.extend(result)
        else:
            # 单线程处理
            for xic_entry in tqdm(xic_entries, total=len(xic_entries), desc="提取XIC"):
                xics.append(self._batch_extract_xic([xic_entry])[0])
        
        return xics

    def _format_xic_entries(self, df: pd.DataFrame) -> list[dict]:
        """格式化xic_entries, 存储提取xic必要的信息"""

        def calculate_isotope_mzs(mz: float, charge: int, num_isotopes: int = 4) -> List[float]:
            """计算同位素峰的质荷比"""
            isotope_mzs = [mz]
            for i in range(1, num_isotopes):
                # 碳同位素间隔约为 1.0033 Da
                isotope_mz = mz + (i * 1.00335483507) / charge
                isotope_mzs.append(isotope_mz)
            return isotope_mzs

        xic_entries = []
        for index, row in tqdm(df.iterrows(), total=len(df), desc="构建XIC Data"):
            precursor_mzs = calculate_isotope_mzs(float(row["precursor_mz"]), int(row["charge"]))
            fragment_mzs = [float(mz) for mz in row["fragment_mz"].split(",")]

            if self.mode == 'rt_range':
                assert 'rt_start' in row and 'rt_stop' in row, "rt_range mode requires rt_start and rt_stop columns"
                rt_start = float(row["rt_start"]) * 60
                rt_stop = float(row["rt_stop"]) * 60
                ms_cluster = self._get_ms_scan_by_rt_range((rt_start, rt_stop), precursor_mzs[0])
            elif self.mode == 'scan_window':
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

    def _batch_extract_xic(self, xic_entries: list[dict]) -> list[tuple[List[XICResult], List[XICResult]]]: 
        """从mz簇中提取XIC"""
        def extract_xic_from_peaks(peaks: list[list[tuple]], mz: float, ppm_tolerance: float) -> tuple[List[float], List[float]]:
            """从峰列表中提取XIC"""
            def binary_search(peaks: list[tuple], mz: float, ppm_tolerance: float) -> tuple[float, float]:
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

            def linear_interpolation(list: list[float]) -> list[float]:
                """线性插值处理缺失值（intensity为None的点）"""
                for i in range(len(list)):
                    if list[i] is None:
                        left_idx = i - 1
                        right_idx = i + 1
                        while left_idx >= 0 and list[left_idx] is None:
                            left_idx -= 1
                        while right_idx < len(list) and list[right_idx] is None:
                            right_idx += 1
                        if left_idx < 0 and right_idx >= len(list):
                            raise ValueError("无法插值：所有数据点都是None")
                        elif left_idx < 0 or right_idx >= len(list):
                            if left_idx >= 0:
                                list[i] = list[left_idx]
                            else:
                                list[i] = list[right_idx]
                        else:
                            list[i] = list[left_idx] + (list[right_idx] - list[left_idx]) * (i - left_idx) / (right_idx - left_idx)
                return list

            intensity_list = []
            ppm_list = []
            
            for peak in peaks:
                if not peak:
                    intensity_list.append(None)
                    ppm_list.append(None)
                else:
                    intensity, ppm = binary_search(peak, mz, ppm_tolerance)
                    intensity_list.append(intensity)
                    ppm_list.append(ppm)
            
            # 线性插值处理缺失值（intensity为None的点）
            if None in intensity_list:
                intensity_list = linear_interpolation(intensity_list)
                ppm_list = linear_interpolation(ppm_list)

            return intensity_list, ppm_list

        xics = []
        for xic_entry in xic_entries:
            precursor_mzs = xic_entry["precursor_mzs"]
            fragment_mzs = xic_entry["fragment_mzs"]
            ms_cluster = xic_entry["ms_cluster"]
            rt_list = [cluster['rt'] for cluster in ms_cluster]

            #precursor_mz的XIC
            precursor_xics = []
            ms1_specs = [cluster['ms1'] for cluster in ms_cluster]
            for precursor_mz in precursor_mzs:
                precursor_intensity_list, precursor_ppm_list = extract_xic_from_peaks(ms1_specs, precursor_mz, self.ppm_tolerance)
                precursor_xics.append(XICResult(rt_list, precursor_intensity_list, precursor_mz, precursor_ppm_list))

            # 获取所有对应范围的ms2谱图
            ms2_specs = [cluster['ms2'] for cluster in ms_cluster]
            #fragment_mz的XIC
            fragment_xics = []
            for fragment_mz in fragment_mzs:
                if all(ms2_spec is None for ms2_spec in ms2_specs):
                    fragment_xics.append(XICResult(rt_list, [0] * len(rt_list), fragment_mz, [self.ppm_tolerance] * len(rt_list)))
                else:
                    fragment_intensity_list, fragment_ppm_list = extract_xic_from_peaks(ms2_specs, fragment_mz, self.ppm_tolerance)
                    fragment_xics.append(XICResult(rt_list, fragment_intensity_list, fragment_mz, fragment_ppm_list))

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
        
        # 如果有足够的XIC，计算平均边界
        if len(all_boundaries) > 0:
            avg_start = int(sum([b[0] for b in all_boundaries]) / len(all_boundaries))
            avg_stop = int(sum([b[1] for b in all_boundaries]) / len(all_boundaries))
        else:
            avg_start = len(xic.rt_array)//2 - self.min_scans//2
            avg_stop = len(xic.rt_array)//2 + self.min_scans//2
        
        # 调整所有XIC
        adjusted_precursor_xics = []
        for xic in precursor_xics:
            if len(xic.rt_array) > 0:
                start_idx = max(0, avg_start)
                stop_idx = min(len(xic.rt_array) - 1, avg_stop)
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
                start_idx = max(0, avg_start)
                stop_idx = min(len(xic.rt_array) - 1, avg_stop)
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
        while bin_index not in self.rt_indices:
            bin_index -= 1
            if bin_index < min(self.rt_indices.keys()):
                break
        if bin_index < min(self.rt_indices.keys()):
            rt_low_index = self.rt_indices[min(self.rt_indices.keys())][0]
        else:
            rt_low_index = self.rt_indices[bin_index][0]

        bin_index = rt_bin
        while bin_index not in self.rt_indices:
            bin_index += 1
            if bin_index > max(self.rt_indices.keys()):
                break
        if bin_index > max(self.rt_indices.keys()):
            rt_high_index = self.rt_indices[max(self.rt_indices.keys())][1]
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
                filtered_ms2 = cluster['ms2'][greedy_index]['peaks']
            else:
                for index, ms2 in enumerate(cluster['ms2']):
                    if ms2['mz_min'] <= precursor_mz <= ms2['mz_max']:
                        filtered_ms2 = ms2['peaks']
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
            if bin_index in self.rt_indices:
                if low_index == -1:
                    low_index = self.rt_indices[bin_index][0]
                    break
        
        for bin_index in range(bin_high, bin_low - 1, -1):
            if bin_index in self.rt_indices:
                if high_index == -1:
                    high_index = self.rt_indices[bin_index][1]
                    break
        
        if low_index == -1 or high_index == -1:
            results = []
        else:
            clusters = self.ms_clusters[low_index:high_index + 1]
            greedy_index = 0
            results = []
            for cluster in clusters:
                result = {'rt': cluster['rt'], 'ms1': cluster['ms1'], 'ms2': None}
                filtered_ms2 = None
                if 0 <= greedy_index < len(cluster['ms2']) and cluster['ms2'][greedy_index]['mz_min'] <= precursor_mz <= cluster['ms2'][greedy_index]['mz_max']:
                    filtered_ms2 = cluster['ms2'][greedy_index]['peaks']
                else:
                    for index, ms2 in enumerate(cluster['ms2']):
                        if ms2['mz_min'] <= precursor_mz <= ms2['mz_max']:
                            filtered_ms2 = ms2['peaks']
                            greedy_index = index
                            break
                result['ms2'] = filtered_ms2
                results.append(result)
        return results
    
    def _format_ms_clusters(self, ms_objects: list[MSObject]):
        """
        将 MSObject 列表转换为簇格式, 并记录每个cluster的rt
        
        参数:
            ms_objects: MSObject 列表
        """
        # 按scan number排序
        ms_objects.sort(key=lambda x: x.retention_time)

        current_cluster = None
        for ms_obj in ms_objects:
            if ms_obj.level == 1:
                if not current_cluster:
                    current_cluster = {'rt': ms_obj.retention_time, 'ms1': ms_obj.peaks, 'ms2': []}
                else:
                    self.ms_clusters.append(current_cluster)
                    rt_index = int(current_cluster['rt'] / self.rt_bin_size)
                    if rt_index not in self.rt_indices:
                        self.rt_indices[rt_index] = [len(self.ms_clusters) - 1, len(self.ms_clusters) - 1]
                    self.rt_indices[rt_index][1] = len(self.ms_clusters) - 1
                    current_cluster = {'rt': ms_obj.retention_time, 'ms1': ms_obj.peaks, 'ms2': []}
            else:
                current_cluster['ms2'].append({'mz': ms_obj.precursor.mz, 'mz_min': ms_obj.precursor.isolation_window[0], 'mz_max': ms_obj.precursor.isolation_window[1], 'peaks': ms_obj.peaks})
        
        if current_cluster:
            self.ms_clusters.append(current_cluster)
        
    