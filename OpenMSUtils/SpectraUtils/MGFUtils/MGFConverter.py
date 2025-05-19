from typing import Any
from ..MSObject import MSObject
from .MGFObject import MGFSpectrum

class MGFConverter:
    @staticmethod
    def to_msobject(spectrum: MGFSpectrum) -> MSObject:
        """
        将MGF的Spectrum对象转换为MSObject
        
        Args:
            spectrum: MGFSpectrum对象
            
        Returns:
            MSObject对象
        """
        # 创建MSObject
        ms_object = MSObject()
        
        # 设置MS级别（MGF通常是MS2）
        ms_object.set_level(2)
        
        # 设置scan信息
        # MGF没有明确的scan number，使用0作为默认值
        scan_number = 0
        
        # 从title中尝试提取scan number
        if spectrum.title and "scan=" in spectrum.title.lower():
            try:
                scan_parts = spectrum.title.lower().split("scan=")[1].split()[0]
                scan_number = int(scan_parts)
            except (ValueError, IndexError):
                pass
        
        ms_object.set_scan(scan_number=scan_number, retention_time=spectrum.rtinseconds)
        
        # 设置前体离子信息
        ms_object.set_precursor(mz=spectrum.pepmass, charge=spectrum.charge)
        
        # 添加峰值
        for mz, intensity in spectrum.peaks:
            ms_object.add_peak(mz, intensity)
        ms_object.sort_peaks() 

        # 添加额外信息
        for key, value in spectrum.additional_info.items():
            ms_object.set_additional_info(key, value)
        
        # 如果有title，添加为额外信息
        if spectrum.title:
            ms_object.set_additional_info("TITLE", spectrum.title)
        
        return ms_object
    
    @staticmethod
    def from_msobject(ms_object: MSObject) -> MGFSpectrum:
        """
        将MSObject转换为MGF的Spectrum对象
        
        Args:
            ms_object: MSObject对象
            
        Returns:
            MGFSpectrum对象
        """
        # 创建MGFSpectrum
        mgf_spectrum = MGFSpectrum()
        
        # 设置标题（使用scan number）
        mgf_spectrum.title = f"Scan={ms_object.scan.scan_number}"
        
        # 设置肽质量（前体离子m/z）
        if ms_object.precursor and ms_object.precursor.mz > 0:
            mgf_spectrum.pepmass = ms_object.precursor.mz
        
        # 设置电荷
        if ms_object.precursor and ms_object.precursor.charge != 0:
            mgf_spectrum.charge = ms_object.precursor.charge
        
        # 设置保留时间
        if ms_object.scan and ms_object.scan.retention_time > 0:
            mgf_spectrum.rtinseconds = ms_object.scan.retention_time
        
        # 添加峰值
        for mz, intensity in ms_object.peaks:
            mgf_spectrum.add_peak(mz, intensity)
        
        # 添加额外信息
        for key, value in ms_object.additional_info.items():
            if key != "TITLE":  # 避免重复添加标题
                mgf_spectrum.set_additional_info(key, value)
        
        return mgf_spectrum 