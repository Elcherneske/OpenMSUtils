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
    rt_array: np.ndarray  # 保留时间数组
    intensity_array: np.ndarray  # 强度数组
    mz: float  # 目标质荷比
    ppm_error: float  # PPM 误差

class XICSExtractor:
    """XIC 提取器类"""
    def __init__(self, mzml_file: str, ppm_tolerance: float = 10.0, rt_bin_size: float = 1.0, num_threads: int = 1):
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
        xic_entries = []
        for index, row in tqdm(df.iterrows(), total=len(df), desc="构建XIC Data"):
            precursor_mzs = self._calculate_isotope_mzs(float(row["precursor_mz"]), int(row["charge"]))
            rt_start = float(row["rt_start"]) * 60
            rt_stop = float(row["rt_stop"]) * 60
            rt = float(row["rt"]) * 60
            fragment_mzs = [float(mz) for mz in row["fragment_mz"].split(",")]
            ms_cluster = self._get_cluster_by_rt_precursor((rt_start, rt_stop), precursor_mzs[0])
            xic_entries.append(
                {
                    "precursor_mzs": precursor_mzs,
                    "fragment_mzs": fragment_mzs,
                    "ms_cluster": ms_cluster
                }
            )

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
                    pool.imap(self._extract_xic_from_mz_cluster, chunks),
                    total=len(chunks),
                    desc=f"提取XIC (使用{self.num_threads}个进程)"
                ))
                
                # 合并结果
                for result in results:
                    xics.extend(result)
        else:
            # 单线程处理
            for xic_entry in tqdm(xic_entries, total=len(xic_entries), desc="提取XIC"):
                xics.append(self._extract_xic_from_mz_cluster([xic_entry])[0])
        
        return xics

    def _extract_xic_from_mz_cluster(self, xic_entries: list[dict]) -> list[tuple[List[XICResult], List[XICResult]]]:
        """从mz簇中提取XIC"""

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
                precursor_intensity_list, precursor_ave_ppm = self._extract_xic_from_peaks(ms1_specs, precursor_mz, self.ppm_tolerance)
                precursor_xics.append(XICResult(np.array(rt_list), np.array(precursor_intensity_list), precursor_mz, precursor_ave_ppm))

            # 获取所有对应范围的ms2谱图
            ms2_specs = [cluster['ms2'] for cluster in ms_cluster]
            #fragment_mz的XIC
            fragment_xics = []
            for fragment_mz in fragment_mzs:
                if all(ms2_spec is None for ms2_spec in ms2_specs):
                    fragment_xics.append(XICResult(np.array(rt_list), np.array([0] * len(rt_list)), fragment_mz, 0))
                else:
                    fragment_intensity_list, fragment_ave_ppm = self._extract_xic_from_peaks(ms2_specs, fragment_mz, self.ppm_tolerance)
                    fragment_xics.append(XICResult(np.array(rt_list), np.array(fragment_intensity_list), fragment_mz, fragment_ave_ppm))

            xics.append((precursor_xics, fragment_xics))
        return xics

    
    def _extract_xic_from_peaks(self, peaks: list[list[tuple]], mz: float, ppm_tolerance: float) -> tuple[List[float], float]:
        """从峰列表中提取XIC"""
        intensity_list = []
        sum_ppm = 0
        valid_ppm_count = 0
        
        for peak in peaks:
            if not peak:
                intensity_list.append(None)
            else:
                intensity, ppm = self._binary_search(peak, mz, ppm_tolerance)
                intensity_list.append(intensity)
                if ppm > 0:  # 只有找到有效峰值时才累加ppm
                    sum_ppm += ppm
                    valid_ppm_count += 1
        
        # 计算平均ppm
        if valid_ppm_count > 0:
            ave_ppm = sum_ppm / valid_ppm_count
        else:
            ave_ppm = 0
        
        # 线性插值处理缺失值（intensity为None的点）
        if None in intensity_list:
            for i in range(len(intensity_list)):
                if intensity_list[i] is None:
                    left_idx = i - 1
                    right_idx = i + 1
                    
                    # 查找左侧非None值
                    while left_idx >= 0 and intensity_list[left_idx] is None:
                        left_idx -= 1
                    
                    # 查找右侧非None值
                    while right_idx < len(intensity_list) and intensity_list[right_idx] is None:
                        right_idx += 1
                    
                    # 如果左右两侧都找不到有效值，则无法插值
                    if left_idx < 0 and right_idx >= len(intensity_list):
                        print(f"mz: {mz}")
                        raise ValueError("无法插值：所有数据点都是None")
                    elif left_idx < 0 or right_idx >= len(intensity_list):
                        # 如果只有一侧有值，使用该侧的值
                        if left_idx >= 0:
                            intensity_list[i] = intensity_list[left_idx]
                        else:
                            intensity_list[i] = intensity_list[right_idx]
                    else:
                        # 线性插值
                        left_val = intensity_list[left_idx]
                        right_val = intensity_list[right_idx]
                        intensity_list[i] = left_val + (right_val - left_val) * (i - left_idx) / (right_idx - left_idx)

        return intensity_list, ave_ppm
    
    def _binary_search(self, peaks: list[tuple], mz: float, ppm_tolerance: float) -> tuple[float, float]:
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
            return 0, 0
    
    def _get_cluster_by_rt_precursor(self, rt_range: tuple[float, float], precursor_mz: float) -> list[dict]:
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
                if greedy_index >= 0 and greedy_index < len(cluster['ms2']) and cluster['ms2'][greedy_index]['mz_min'] <= precursor_mz <= cluster['ms2'][greedy_index]['mz_max']:
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
        将 MSObject 列表转换为簇格式
        
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
        
    def _calculate_isotope_mzs(self, mz: float, charge: int, num_isotopes: int = 4) -> List[float]:
        """计算同位素峰的质荷比"""
        isotope_mzs = [mz]
        for i in range(1, num_isotopes):
            # 碳同位素间隔约为 1.0033 Da
            isotope_mz = mz + (i * 1.00335483507) / charge
            isotope_mzs.append(isotope_mz)
        return isotope_mzs
    