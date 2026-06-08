from MSUtils import SpectraObject
from typing import Dict, List, Tuple
import numpy as np

class IonMobilityUtils:
    def __init__(self):
        super().__init__()
    
    @staticmethod
    def _merge_peaks_by_mz(peaks: np.ndarray, mz_tolerance: float) -> List[Tuple[float, float]]:
        """
        根据m/z容差合并相似m/z值的峰强度
        
        参数:
            peaks: 峰列表，每个峰为np.ndarray
            mz_tolerance: m/z容差，单位为ppm
            
        返回:
            合并后的峰列表
        """
        if not peaks or len(peaks) == 0:
            return np.array([])
        
        peaks = peaks[np.argsort(peaks[:, 0])]
        merged_peaks = []
        current_group = [peaks[0]]
        current_mz = peaks[0][0]
        
        for peak in peaks[1:]:
            mz = peak[0]
            intensity = peak[1]
            # 计算当前峰与当前组m/z的差异（ppm）
            tolerance_da = current_mz * mz_tolerance / 1e6
            
            if abs(mz - current_mz) <= tolerance_da:
                current_group.append(peak)
            else:
                if current_group:
                    total_intensity = sum(p[1] for p in current_group)
                    weighted_mz = sum(p[0] * p[1] for p in current_group) / total_intensity if total_intensity > 0 else current_mz
                    merged_peaks.append((weighted_mz, total_intensity))
                
                # 开始新组
                current_group = [peak]
                current_mz = mz
        
        # 处理最后一组
        if current_group:
            total_intensity = sum(p[1] for p in current_group)
            weighted_mz = sum(p[0] * p[1] for p in current_group) / total_intensity if total_intensity > 0 else current_mz
            merged_peaks.append((weighted_mz, total_intensity))
        
        return np.array(merged_peaks, dtype=np.float32)
    
    @staticmethod
    def parse_ion_mobility(spectra_objects: List[SpectraObject], rt_range=None, mz_tolerance=10, rt_tolerance=None) -> Dict[float, list[tuple[float, float]]]:  
        """
        解析淌度谱数据并返回字典
        
        参数:
            spectra_objects: SpectraObject对象列表
            rt_range: 保留时间范围，格式为(min_rt, max_rt)
            mz_tolerance: m/z容差，默认为10ppm
            rt_tolerance: 保留时间容差，默认为None
            
        返回:
            ion_mobility_spectrum: 淌度谱字典，键为漂移时间，值为对应的峰列表
        """
        # 如果提供了rt_range，筛选在保留时间范围内的MSObject
        filtered_spectra_objects = spectra_objects
        if rt_range is not None:
            min_rt, max_rt = rt_range
            filtered_spectra_objects = [spectra_obj for spectra_obj in spectra_objects if spectra_obj.retention_time is not None and min_rt <= spectra_obj.retention_time <= max_rt]
        
        # 初始化淌度谱字典
        ion_mobility_spectrum = {}
        # 解析淌度谱数据
        for spectra_object in filtered_spectra_objects:
            drift_time = spectra_object.drift_time
            if not drift_time or drift_time < 0:
                continue

            # 如果提供了rt_tolerance，查找在容差范围内的已有drift_time
            found_key = None
            if rt_tolerance is not None:
                for existing_dt in ion_mobility_spectrum.keys():
                    if abs(existing_dt - drift_time) <= rt_tolerance:
                        found_key = existing_dt
                        break
            
            # 确定使用的键
            key = found_key if found_key is not None else drift_time
            
            # 初始化字典项（如果不存在）
            if key not in ion_mobility_spectrum:
                ion_mobility_spectrum[key] = []
                
            ion_mobility_spectrum[key].extend(spectra_object.peaks)
        
        # 对每个漂移时间的峰列表根据mz_tolerance合并相似m/z值的峰
        for drift_time, peaks in ion_mobility_spectrum.items():
            ion_mobility_spectrum[drift_time] = IonMobilityUtils._merge_peaks_by_mz(peaks, mz_tolerance)
        
        return ion_mobility_spectrum


