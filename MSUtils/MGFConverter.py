import re
from .MSObject import SpectraObject

class MGFConverter:
    @staticmethod
    def to_spectra_object(lines: list[str]) -> SpectraObject:
        """
        将MS1/MS2的Spectrum对象转换为SpectraObject
        
        Args:
            lines: 包含MS1/MS2数据的行列表
            
        Returns:    
            SpectraObject对象
        """
        # 创建SpectraObject
        spectra_object = SpectraObject()

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
            if not line or line.startswith('#') or line.startswith('BEGIN IONS') or line.startswith('END IONS'):
                continue

            if line.startswith('TITLE'):
                if "scan=" in line:
                    scan_number_str = line.split("scan=")[-1]
                    match = re.search(r'(\d+)', scan_number_str)
                    if match:
                        scan_number = int(match.group(1))
                continue
            
            # 开始新的谱图
            if line.startswith('RTINSECONDS'):
                parts = line.split('=')
                if not len(parts) == 2:
                    continue
                retention_time = float(parts[1])
                continue
            
            # 处理信息行
            if line.startswith('PEPMASS'):
                parts = line.split('=')
                if not len(parts) == 2:
                    continue
                precursor_mz = float(parts[1].split()[0])
                continue
            
            # 处理电荷行（仅MS2）
            if line.startswith('CHARGE'):
                parts = line.split('=')
                if not len(parts) == 2:
                    continue
                if parts[1].endswith('+'):
                    precursor_charge = int(parts[1][:-1])
                elif parts[1].endswith('-'):
                    precursor_charge = -int(parts[1][:-1])
                else:
                    precursor_charge = int(parts[1])
                continue
            
            # 处理峰值数据
            parts = line.split()
            if not len(parts) >= 2:
                continue
            mz = float(parts[0])
            intensity = float(parts[1])
            peaks.append((mz, intensity))
        
        spectra_object.set_level(ms_level)
        spectra_object.set_scan(scan_number=scan_number, retention_time=retention_time, drift_time=drift_time, scan_window=scan_window)
        spectra_object.set_precursor(mz=precursor_mz, charge=precursor_charge, ref_scan_number=-1, activation_method=activation_method, activation_energy=activation_energy, isolation_window=isolation_window)
        spectra_object.set_peaks(peaks)

        return spectra_object
    
    @staticmethod
    def from_spectra_object(spectra_object: SpectraObject) -> list[str]:
        """
        将SpectraObject转换为MGF的Spectrum对象
        
        Args:
            spectra_object: SpectraObject对象
            
        Returns:
            包含MGF数据的行列表
        """
        lines = []
        lines.append(f"BEGIN IONS")
        lines.append(f"TITLE=Scan={spectra_object.scan_number}")
        lines.append(f"RTINSECONDS={spectra_object.retention_time}")
        lines.append(f"PEPMASS={spectra_object.precursor_mz}")
        lines.append(f"CHARGE={spectra_object.precursor_charge}")
        
        for mz, intensity in spectra_object.peaks:
            lines.append(f"{mz}\t{intensity}")
        
        lines.append(f"END IONS")
        
        return lines