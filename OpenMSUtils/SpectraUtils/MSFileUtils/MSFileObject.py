class MSSpectrum(object):
    """
    MS1/MS2格式的质谱数据对象
    """
    def __init__(self, level=1):
        self._level = level  # MS级别，1表示MS1，2表示MS2
        self._scan_number = 0
        self._retention_time = 0.0
        self._precursor_mz = 0.0
        self._precursor_charge = 0
        self._peaks = []
        self._additional_info = {}
    
    @property
    def level(self):
        """获取MS级别"""
        return self._level
    
    @level.setter
    def level(self, value):
        """设置MS级别"""
        self._level = value
    
    @property
    def scan_number(self):
        """获取扫描编号"""
        return self._scan_number
    
    @scan_number.setter
    def scan_number(self, value):
        """设置扫描编号"""
        self._scan_number = value
    
    @property
    def retention_time(self):
        """获取保留时间"""
        return self._retention_time
    
    @retention_time.setter
    def retention_time(self, value):
        """设置保留时间"""
        self._retention_time = value
    
    @property
    def precursor_mz(self):
        """获取前体离子m/z"""
        return self._precursor_mz
    
    @precursor_mz.setter
    def precursor_mz(self, value):
        """设置前体离子m/z"""
        self._precursor_mz = value
    
    @property
    def precursor_charge(self):
        """获取前体离子电荷"""
        return self._precursor_charge
    
    @precursor_charge.setter
    def precursor_charge(self, value):
        """设置前体离子电荷"""
        self._precursor_charge = value
    
    @property
    def peaks(self):
        """获取峰值列表"""
        return self._peaks
    
    @peaks.setter
    def peaks(self, value):
        """设置峰值列表"""
        self._peaks = value
    
    @property
    def additional_info(self):
        """获取额外信息字典"""
        return self._additional_info
    
    @additional_info.setter
    def additional_info(self, value):
        """设置额外信息字典"""
        self._additional_info = value
    
    def add_peak(self, mz, intensity):
        """添加峰值"""
        self._peaks.append((mz, intensity))
    
    def set_additional_info(self, key, value):
        """设置额外信息"""
        self._additional_info[key] = value
    
    def to_ms_string(self):
        """转换为MS1/MS2格式字符串"""
        lines = []
        
        # 添加扫描信息行
        if self._level == 1:
            lines.append(f"S\t{self._level}\t{self._scan_number}")
        else:
            lines.append(f"S\t{self._level}\t{self._scan_number}\t{self._precursor_mz}")
        
        # 添加保留时间
        if self._retention_time > 0:
            lines.append(f"I\tRTime\t{self._retention_time}")
        
        # 添加额外信息
        for key, value in self._additional_info.items():
            lines.append(f"I\t{key}\t{value}")
        
        # 添加电荷信息（仅MS2）
        if self._level == 2 and self._precursor_charge != 0:
            lines.append(f"Z\t{self._precursor_charge}\t{self._precursor_mz * abs(self._precursor_charge)}")
        
        # 添加峰值数据
        for mz, intensity in self._peaks:
            lines.append(f"{mz} {intensity}")
        
        return "\n".join(lines)

class MSFileObject(object):
    """
    MS1/MS2文件对象，包含多个MSSpectrum
    """
    def __init__(self, level=1):
        self._level = level  # 文件级别，1表示MS1，2表示MS2
        self._spectra = []
        self._metadata = {}
    
    @property
    def level(self):
        """获取文件级别"""
        return self._level
    
    @level.setter
    def level(self, value):
        """设置文件级别"""
        self._level = value
    
    @property
    def spectra(self):
        """获取谱图列表"""
        return self._spectra
    
    @spectra.setter
    def spectra(self, value):
        """设置谱图列表"""
        self._spectra = value
    
    @property
    def metadata(self):
        """获取元数据字典"""
        return self._metadata
    
    @metadata.setter
    def metadata(self, value):
        """设置元数据字典"""
        self._metadata = value
    
    def add_spectrum(self, spectrum):
        """添加谱图"""
        self._spectra.append(spectrum)
    
    def set_metadata(self, key, value):
        """设置元数据"""
        self._metadata[key] = value
    
    def to_ms_string(self):
        """转换为MS1/MS2格式字符串"""
        lines = []
        
        # 添加元数据作为头信息
        for key, value in self._metadata.items():
            lines.append(f"H\t{key}\t{value}")
        
        # 添加空行分隔头信息和谱图数据
        lines.append("")
        
        # 添加所有谱图
        for spectrum in self._spectra:
            lines.append(spectrum.to_ms_string())
            lines.append("")  # 添加空行分隔谱图
        
        return "\n".join(lines) 

