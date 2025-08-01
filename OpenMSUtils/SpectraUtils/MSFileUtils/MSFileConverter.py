import re
from typing import Any
from ..MSObject import MSObject
from .MSFileObject import MSSpectrum

class MSFileConverter:
    @staticmethod
    def to_msobject(lines: list[str]) -> MSObject:
        """
        将MS1/MS2的Spectrum对象转换为MSObject
        
        Args:
            lines: 包含MS1/MS2数据的行列表
            
        Returns:    
            MSObject对象
        """
        # 创建MSObject
        ms_object = MSObject()

        ms_level = 1
        scan_number = -1
        retention_time = 0.0
        drift_time = 0.0
        scan_window = (0.0, 0.0)
        precursor_mz = 0.0
        precursor_charge = 0
        activation_method = 'unknown'
        activation_energy = 0.0
        isolation_window = (0.0, 0.0)
        peaks = []

        for line in lines:
            line = line.strip()
            
            # 跳过空行
            if not line or line.startswith('#') or line.startswith('H'):
                continue
            
            # 开始新的谱图
            if line.startswith('S'):
                parts = line.split()
                if len(parts) == 4:
                    ms_level = 2
                    scan_number = int(parts[1])
                    precursor_mz = float(parts[3])
                elif len(parts) == 3:
                    ms_level = 1
                    scan_number = int(parts[1])
                continue
            
            # 处理信息行
            if line.startswith('I'):
                parts = line.split()
                if len(parts) >= 3:
                    if parts[1] == "RTime":
                        retention_time = float(parts[2])
                continue
            
            # 处理电荷行（仅MS2）
            if line.startswith('Z') and ms_level == 2:
                parts = line.split()
                if len(parts) >= 2:
                    precursor_charge = int(parts[1])
                continue
            
            # 处理峰值数据
            parts = line.split()
            if not len(parts) >= 2:
                continue
            mz = float(parts[0])
            intensity = float(parts[1])
            peaks.append((mz, intensity))
        
        ms_object.set_level(ms_level)
        ms_object.set_scan(scan_number=scan_number, retention_time=retention_time, drift_time=drift_time, scan_window=scan_window)
        ms_object.set_precursor(mz=precursor_mz, charge=precursor_charge, ref_scan_number=-1, activation_method=activation_method, activation_energy=activation_energy, isolation_window=isolation_window)
        ms_object.set_peaks(peaks)

        return ms_object
    
    @staticmethod
    def from_msobject(ms_object: MSObject) -> list[str]:
        """
        将MSObject转换为MS1/MS2的Spectrum对象
        
        Args:
            ms_object: MSObject对象
            
        Returns:
            包含MS1/MS2数据的行列表
        """
        lines = []
        if ms_object.level == 1:
            lines.append(f"S\t{ms_object.scan_number}\t{ms_object.scan_number}")
        elif ms_object.level == 2:
            lines.append(f"S\t{ms_object.scan_number}\t{ms_object.scan_number}\t{ms_object.precursor_mz}")
        lines.append(f"I\tRTime\t{ms_object.retention_time}")

        if ms_object.level == 2:
            lines.append(f"Z\t{ms_object.precursor_charge}")
        
        for mz, intensity in ms_object.peaks:
            lines.append(f"{mz}\t{intensity}")

        return lines