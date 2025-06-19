import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from .MSObject import MSObject
from .XICSExtractor import XICResult

class SpectraPlotter:
    def __init__(self, mode='plt'):
        self.mode = mode
        pass

    def plot_mz(self, ms_object: MSObject):
        """
        绘制 m/z 图
        
        参数:
            ms_object: MSObject 对象
        
        返回:
            None
        """

        def plot_plt(x, y, title, xlabel, ylabel, color='blue'):
            plt.figure(figsize=(10, 6))
            plt.stem(x, y, markerfmt=" ", linefmt=color, basefmt=color)
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.show()

        def plot_plotly(x, y, title, xlabel, ylabel, color='blue'):
            fig = go.Figure()
            
            # 添加针状线条
            for x_val, y_val in zip(x, y):
                fig.add_trace(go.Scatter(
                    x=[x_val, x_val],
                    y=[0, y_val],
                    mode='lines',
                    line=dict(color=color, width=1),
                    showlegend=False,
                    hoverinfo='skip'
                ))
            
            # 添加数据点
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode='markers',
                marker=dict(
                    size=1,
                    color=color,
                    symbol='circle'
                ),
                name='peaks',
                hovertemplate='<b>m/z:</b> %{x}<br><b>Intensity:</b> %{y}<extra></extra>'
            ))
            
            fig.update_layout(
                title=title, 
                xaxis_title=xlabel, 
                yaxis_title=ylabel,
                showlegend=True,
                hovermode='closest'
            )

            fig.show()
            return fig

        peaks = ms_object.peaks
        if len(peaks) == 0:
            print("No peaks found in the spectrum")
            return
        # 提取 mz 和 intensity 数据
        mz_list = [peak[0] for peak in peaks]
        intensity_list = [peak[1] for peak in peaks]
        max_intensity = max(intensity_list) if intensity_list else 1
        intensity_list = [intensity / max_intensity for intensity in intensity_list]

        if self.mode == 'plt':
            plot_plt(mz_list, intensity_list, "Mass Spectrum", "m/z", "Intensity")
        elif self.mode == 'plotly':
            plot_plotly(mz_list, intensity_list, "Mass Spectrum", "m/z", "Intensity")

    def plot_mz_comparison(self, ms_object1: MSObject, ms_object2: MSObject):
        """
        绘制两个 m/z 图的对比图
        
        参数:
            ms_object1: 第一个 MSObject 对象
            ms_object2: 第二个 MSObject 对象
        
        返回:
            None
        """

        def plot_plt_comparison(x1, y1, x2, y2, xlabel, ylabel):
            fig, ax = plt.subplots(1, 1, figsize=(12, 8))
            
            # 绘制第一个谱图（正半轴）
            ax.stem(x1, y1, markerfmt=" ", linefmt='blue', basefmt='blue', label="Spectrum 1")
            
            # 绘制第二个谱图（负半轴）
            ax.stem(x2, [-y for y in y2], markerfmt=" ", linefmt='red', basefmt='red', label="Spectrum 2")
            
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            
            # 设置y轴刻度为正值
            y_ticks = ax.get_yticks()
            ax.set_yticklabels([abs(tick) for tick in y_ticks])
            
            plt.tight_layout()
            plt.show()

        def plot_plotly_comparison(x1, y1, x2, y2, xlabel, ylabel):
            fig = go.Figure()
            
            # 添加第一个谱图的针状线条（正半轴）
            for x_val, y_val in zip(x1, y1):
                fig.add_trace(go.Scatter(
                    x=[x_val, x_val],
                    y=[0, y_val],
                    mode='lines',
                    line=dict(color='blue', width=1),
                    showlegend=False,
                    hoverinfo='skip'
                ))
            
            # 添加第一个谱图的数据点
            fig.add_trace(go.Scatter(
                x=x1,
                y=y1,
                mode='markers',
                marker=dict(size=1, color='blue', symbol='circle'),
                name="Spectrum 1",
                hovertemplate='<b>m/z:</b> %{x}<br><b>Intensity:</b> %{y}<extra></extra>'
            ))
            
            # 添加第二个谱图的针状线条（负半轴）
            for x_val, y_val in zip(x2, y2):
                fig.add_trace(go.Scatter(
                    x=[x_val, x_val],
                    y=[0, -y_val],
                    mode='lines',
                    line=dict(color='red', width=1),
                    showlegend=False,
                    hoverinfo='skip'
                ))
            
            # 添加第二个谱图的数据点
            fig.add_trace(go.Scatter(
                x=x2,
                y=[-y for y in y2],
                mode='markers',
                marker=dict(size=1, color='red', symbol='circle'),
                name="Spectrum 2",
                hovertemplate='<b>m/z:</b> %{x}<br><b>Intensity:</b> %{y}<extra></extra>'
            ))
            
            fig.update_layout(
                title="Spectrum 1 vs Spectrum 2",
                xaxis_title=xlabel,
                yaxis_title=ylabel,
                showlegend=True,
                hovermode='closest',
                height=600
            )
            
            # 设置y轴刻度为正值
            fig.update_yaxes(tickformat='.0f', tickmode='auto')

            fig.show()
            return fig

        # 处理第一个 MSObject
        peaks1 = ms_object1.peaks
        if len(peaks1) == 0:
            print("No peaks found in the first spectrum")
            return
            
        mz_list1 = [peak[0] for peak in peaks1]
        intensity_list1 = [peak[1] for peak in peaks1]
        max_intensity1 = max(intensity_list1) if intensity_list1 else 1
        intensity_list1 = [intensity / max_intensity1 for intensity in intensity_list1]

        # 处理第二个 MSObject
        peaks2 = ms_object2.peaks
        if len(peaks2) == 0:
            print("No peaks found in the second spectrum")
            return
            
        mz_list2 = [peak[0] for peak in peaks2]
        intensity_list2 = [peak[1] for peak in peaks2]
        max_intensity2 = max(intensity_list2) if intensity_list2 else 1
        intensity_list2 = [intensity / max_intensity2 for intensity in intensity_list2]

        if self.mode == 'plt':
            plot_plt_comparison(mz_list1, intensity_list1, mz_list2, intensity_list2, "m/z", "Intensity")
        elif self.mode == 'plotly':
            plot_plotly_comparison(mz_list1, intensity_list1, mz_list2, intensity_list2, "m/z", "Intensity")
    
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
        def plot_plt(precursor_xics, fragment_xics):
            fig = plt.figure(figsize=(15, 10))
            ax1 = fig.add_subplot(2, 2, 1)
            ax2 = fig.add_subplot(2, 2, 2)
            ax3 = fig.add_subplot(2, 2, 3)
            ax4 = fig.add_subplot(2, 2, 4)
            
            # 前体离子强度图
            for xic in precursor_xics:
                label = f"(m/z: {xic.mz:.4f}, ppm: {np.average(np.array(xic.ppm_array)):.2f})"
                ax1.plot(xic.rt_array, xic.intensity_array, label=label)
            ax1.set_title("Precursor XIC - Intensity")
            ax1.set_ylabel("Intensity")
            ax1.legend()
            ax1.grid(True, linestyle='--', alpha=0.7)
            
            # 碎片离子强度图
            for xic in fragment_xics:
                label = f"(m/z: {xic.mz:.4f}, ppm: {np.average(np.array(xic.ppm_array)):.2f})"
                ax2.plot(xic.rt_array, xic.intensity_array, label=label)
            ax2.set_title("Fragment XIC - Intensity")
            ax2.set_ylabel("Intensity")
            ax2.legend()
            ax2.grid(True, linestyle='--', alpha=0.7)
            
            # 前体离子ppm图
            for xic in precursor_xics:
                label = f"m/z: {xic.mz:.4f}"
                ax3.plot(xic.rt_array, xic.ppm_array, label=label)
            ax3.set_title("Precursor XIC - PPM")
            ax3.set_xlabel("Retention Time (min)")
            ax3.set_ylabel("PPM")
            ax3.legend()
            ax3.grid(True, linestyle='--', alpha=0.7)
            
            # 碎片离子ppm图
            for xic in fragment_xics:
                label = f"m/z: {xic.mz:.4f}"
                ax4.plot(xic.rt_array, xic.ppm_array, label=label)
            ax4.set_title("Fragment XIC - PPM")
            ax4.set_xlabel("Retention Time (min)")
            ax4.set_ylabel("PPM")
            ax4.legend()
            ax4.grid(True, linestyle='--', alpha=0.7)

            plt.tight_layout()
            return fig
        
        def plot_plotly(precursor_xics, fragment_xics):
            # 创建2x2子图
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=("Precursor XIC - Intensity", "Fragment XIC - Intensity", 
                              "Precursor XIC - PPM", "Fragment XIC - PPM"),
                specs=[[{"secondary_y": False}, {"secondary_y": False}],
                       [{"secondary_y": False}, {"secondary_y": False}]]
            )
            
            # 前体离子强度图
            for xic in precursor_xics:
                label = f"(m/z: {xic.mz:.4f}, ppm: {np.average(np.array(xic.ppm_array)):.2f})"
                fig.add_trace(
                    go.Scatter(
                        x=xic.rt_array,
                        y=xic.intensity_array,
                        mode='lines',
                        name=label,
                        showlegend=True
                    ),
                    row=1, col=1
                )
            
            # 碎片离子强度图
            for xic in fragment_xics:
                label = f"(m/z: {xic.mz:.4f}, ppm: {np.average(np.array(xic.ppm_array)):.2f})"
                fig.add_trace(
                    go.Scatter(
                        x=xic.rt_array,
                        y=xic.intensity_array,
                        mode='lines',
                        name=label,
                        showlegend=True
                    ),
                    row=1, col=2
                )
            
            # 前体离子ppm图
            for xic in precursor_xics:
                label = f"m/z: {xic.mz:.4f}"
                fig.add_trace(
                    go.Scatter(
                        x=xic.rt_array,
                        y=xic.ppm_array,
                        mode='lines',
                        name=label,
                        showlegend=True
                    ),
                    row=2, col=1
                )
            
            # 碎片离子ppm图
            for xic in fragment_xics:
                label = f"m/z: {xic.mz:.4f}"
                fig.add_trace(
                    go.Scatter(
                        x=xic.rt_array,
                        y=xic.ppm_array,
                        mode='lines',
                        name=label,
                        showlegend=True
                    ),
                    row=2, col=2
                )
            
            # 更新布局
            fig.update_layout(
                width=1200,
                height=800,
                title_text="XIC Analysis",
                showlegend=True,
                hovermode='closest'
            )
            
            # 更新坐标轴标签
            fig.update_xaxes(title_text="Retention Time (min)", row=2, col=1)
            fig.update_xaxes(title_text="Retention Time (min)", row=2, col=2)
            fig.update_yaxes(title_text="Intensity", row=1, col=1)
            fig.update_yaxes(title_text="Intensity", row=1, col=2)
            fig.update_yaxes(title_text="PPM", row=2, col=1)
            fig.update_yaxes(title_text="PPM", row=2, col=2)
            
            # 添加网格
            for i in range(1, 3):
                for j in range(1, 3):
                    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', row=i, col=j)
                    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', row=i, col=j)
            
            return fig
        
        if self.mode == 'plt':
            fig = plot_plt(precursor_xics, fragment_xics)
        elif self.mode == 'plotly':
            fig = plot_plotly(precursor_xics, fragment_xics)
        
        if output_file:
            if self.mode == 'plt':
                fig.savefig(output_file, dpi=300)
            elif self.mode == 'plotly':
                fig.write_image(output_file, width=1200, height=800)
            print(f"图表已保存至: {output_file}")
        else:
            fig.show()

if __name__ == "__main__":
    pass



