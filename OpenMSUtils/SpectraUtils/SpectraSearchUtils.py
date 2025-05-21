import numpy as np
from typing import List, Dict, Tuple
from .MSObject import MSObject


class BinnedSpectra:
    def __init__(self, spectra: MSObject|List[Tuple[float, float]], bin_size: float=1.0):
        if isinstance(spectra, MSObject):
            self.spectra = spectra.peaks
        elif isinstance(spectra, list):
            self.spectra = spectra
        else:
            raise TypeError("unsupported spectra type")
        self.spectra.sort(key=lambda x: x[0])
        self.bin_size = bin_size
        self.bin_indices = self._generate_bin_indices(spectra)
    
    def search_peaks(self, mz_range: Tuple[float, float]) -> List[Tuple[float, float]]:
        """
        搜索指定mz范围内的峰值
        
        参数:
            mz_range: 包含mz范围的元组 (mz_low, mz_high)
            
        返回:
            包含mz范围的峰值列表
        """
        mz_low, mz_high = mz_range
        bin_low = int(mz_low / self.bin_size)
        bin_high = int(mz_high / self.bin_size)
        low_index = -1
        high_index = -1
        for bin_index in range(bin_low, bin_high + 1):
            if bin_index in self.bin_indices:
                if low_index == -1:
                    low_index = self.bin_indices[bin_index][0]
                high_index = self.bin_indices[bin_index][1]
        
        if low_index == -1 or high_index == -1:
            result = []
        else:
            result = self.spectra[low_index:high_index + 1]
        result = [peak for peak in result if peak[0] >= mz_low and peak[0] <= mz_high]
        result.sort(key=lambda x: x[0])
        return result

    def _generate_bin_indices(self, data: MSObject|list[Tuple[float, float]]) -> Dict[int, Tuple[int, int]]:
        """
        根据输入的MSObject或者质荷比列表生成bin索引
        
        参数:
            data: MSObject对象或者质荷比列表
            bin_size: bin的大小
            
        返回:
            包含mz到索引映射的字典
        """
        # 提取mz值
        mz_values = []
        if isinstance(data, MSObject):
            if data.peaks is not None:
                mz_values = [peak[0] for peak in data.peaks]
            else:
                raise TypeError("unsupported MSObject type")
        else:
            mz_values = [peak[0] for peak in data]
            
        # 检查mz值是否单调递增
        for i in range(1, len(mz_values)):
            if mz_values[i] < mz_values[i-1]:
                raise ValueError("mz list must be monotonically increasing")
        
        # 生成bin范围和索引
        mz_to_index = {}
        if not mz_values:
            return mz_to_index
        
        for index, mz in enumerate(mz_values):
            bin_index = int(mz / self.bin_size)
            if not bin_index in mz_to_index:
                mz_to_index[bin_index] = (index, index)
            else:
                mz_to_index[bin_index] = (mz_to_index[bin_index][0], index)
            
        return mz_to_index


