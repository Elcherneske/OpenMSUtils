import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional

from .MSObject import MSObject
from .XICSExtractor import XICResult

class SpectraPlotter:
    def __init__(self):
        pass

    def plot_mz(self, ms_object: MSObject):
        peaks = ms_object.peaks
        if len(peaks) == 0:
            print("No peaks found in the spectrum")
            return
        # 提取 mz 和 intensity 数据
        mz_list = [peak[0] for peak in peaks]
        intensity_list = [peak[1] for peak in peaks]
        max_intensity = max(intensity_list) if intensity_list else 1
        intensity_list = [intensity / max_intensity for intensity in intensity_list]

        plt.figure(figsize=(10, 6))
        plt.stem(mz_list, intensity_list, markerfmt=" ")
        # 设置标题和标签
        plt.title("Mass Spectrum")
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.show()
    
    def plot_ion_mobility(self, ion_mobility_spectrum: dict, time_range: tuple[float, float] = None, time_bins: int = 500, mz_range: tuple[float, float] = None, mz_bins: int = 500):
        if time_range is None or mz_range is None:
            min_time = min(ion_mobility_spectrum.keys())
            max_time = max(ion_mobility_spectrum.keys())
            min_mz = min([peak.mz for peaks in ion_mobility_spectrum.values() for peak in peaks])
            max_mz = max([peak.mz for peaks in ion_mobility_spectrum.values() for peak in peaks])
            if time_range is None:
                time_range = (min_time, max_time)
            if mz_range is None:
                mz_range = (min_mz, max_mz)
        time_step = (time_range[1] - time_range[0]) / time_bins
        mz_step = (mz_range[1] - mz_range[0]) / mz_bins
        time_list = np.linspace(time_range[0], time_range[1], time_bins)
        mz_list = np.linspace(mz_range[0], mz_range[1], mz_bins)
        ion_mobility_matrix = np.zeros((len(time_list), len(mz_list)))
        for time, peaks in ion_mobility_spectrum.items():
            for peak in peaks:
                mz = peak.mz
                intensity = peak.intensity
                mz_index = (int)((mz - mz_range[0]) / mz_step)
                time_index = (int)((time - time_range[0]) / time_step)
                if mz_index < 0 or mz_index >= len(mz_list) or time_index < 0 or time_index >= len(time_list):
                    continue
                ion_mobility_matrix[time_index][mz_index] += intensity

        plt.imshow(ion_mobility_matrix, cmap="hot", interpolation='nearest')
        plt.colorbar()
        plt.xlabel("m/z")
        plt.ylabel("Time")
        plt.title("Ion Mobility Spectrum")
        plt.show()
    
    def plot_xics(self, precursor_xics: List[XICResult], fragment_xics: List[XICResult], output_file: Optional[str] = None):
        """绘制 XIC 图
        
        参数:
            precursor_xics: 前体离子 XIC 结果列表
            fragment_xics: 碎片离子 XIC 结果列表
            peptide_info: 肽段信息
            output_file: 输出文件路径，如果为 None 则显示图像
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        # 绘制前体离子 XIC
        for xic in precursor_xics:
            label = f"(m/z: {xic.mz:.4f}, ppm: {np.average(np.array(xic.ppm_array)):.2f})"
            ax1.plot(xic.rt_array, xic.intensity_array, label=label)
            
        ax1.set_title(f"Precursor XIC")
        ax1.set_ylabel("Intensity")
        ax1.legend()
        ax1.grid(True, linestyle='--', alpha=0.7)
        
        # 绘制碎片离子 XIC
        for xic in fragment_xics:
            label = f"(m/z: {xic.mz:.4f}, ppm: {np.average(np.array(xic.ppm_array)):.2f})"
            ax2.plot(xic.rt_array, xic.intensity_array, label=label)
            
        ax2.set_title("Fragment XIC")
        ax2.set_xlabel("Retention Time (min)")
        ax2.set_ylabel("Intensity")
        ax2.legend()
        ax2.grid(True, linestyle='--', alpha=0.7)
        
        # 添加 RT 范围标记
        for ax in [ax1, ax2]:
            ax.axvline(x=precursor_xics[0].rt_array[0], color='r', linestyle='--', alpha=0.5)
            ax.axvline(x=precursor_xics[0].rt_array[-1], color='r', linestyle='--', alpha=0.5)
            ax.axvline(x=precursor_xics[0].rt_array[len(precursor_xics[0].rt_array) // 2], color='g', linestyle='-', alpha=0.5)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300)
            print(f"图表已保存至: {output_file}")
        else:
            plt.show()

if __name__ == "__main__":
    pass



