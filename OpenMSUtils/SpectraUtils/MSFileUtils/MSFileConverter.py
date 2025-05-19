from typing import Any
from ..MSObject import MSObject
from .MSFileObject import MSSpectrum

class MSFileConverter:
    @staticmethod
    def to_msobject(spectrum: MSSpectrum) -> MSObject:
        """
        将MS1/MS2的Spectrum对象转换为MSObject
        
        Args:
            spectrum: MSSpectrum对象
            
        Returns:
            MSObject对象
        """
        # 创建MSObject
        ms_object = MSObject()
        
        # 设置MS级别
        ms_object.set_level(spectrum.level)
        
        # 设置scan信息
        ms_object.set_scan(scan_number=spectrum.scan_number, retention_time=spectrum.retention_time)
        
        # 设置前体离子信息（如果是MS2）
        if spectrum.level == 2:
            ms_object.set_precursor(mz=spectrum.precursor_mz, charge=spectrum.precursor_charge)
        
        # 添加峰值
        for mz, intensity in spectrum.peaks:
            ms_object.add_peak(mz, intensity)
        ms_object.sort_peaks()

        # 添加额外信息
        for key, value in spectrum.additional_info.items():
            ms_object.set_additional_info(key, value)
        
        return ms_object
    
    @staticmethod
    def from_msobject(ms_object: MSObject) -> MSSpectrum:
        """
        将MSObject转换为MS1/MS2的Spectrum对象
        
        Args:
            ms_object: MSObject对象
            
        Returns:
            MSSpectrum对象
        """
        # 创建MSSpectrum
        ms_spectrum = MSSpectrum(level=ms_object.level)
        
        # 设置扫描编号
        ms_spectrum.scan_number = ms_object.scan.scan_number
        
        # 设置保留时间
        ms_spectrum.retention_time = ms_object.scan.retention_time
        
        # 设置前体离子信息（如果是MS2）
        if ms_object.level == 2 and ms_object.precursor:
            ms_spectrum.precursor_mz = ms_object.precursor.mz
            ms_spectrum.precursor_charge = ms_object.precursor.charge
        
        # 添加峰值
        for mz, intensity in ms_object.peaks:
            ms_spectrum.add_peak(mz, intensity)
        
        # 添加额外信息
        for key, value in ms_object.additional_info.items():
            ms_spectrum.set_additional_info(key, value)
        
        return ms_spectrum 