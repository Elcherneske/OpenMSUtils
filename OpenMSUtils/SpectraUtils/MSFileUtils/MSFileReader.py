import os
import re
from tqdm import tqdm
from .MSFileObject import MSFileObject, MSSpectrum
from .MSFileConverter import MSFileConverter

class MSFileReader(object):
    def __init__(self):
        super().__init__()
    
    def read(self, filename):
        """
        读取MS1/MS2文件并解析为MSFileObject对象
        
        Args:
            filename: MS1/MS2文件路径
            
        Returns:
            MSFileObject: 包含MS数据的对象
        """
        if not os.path.exists(filename):
            raise ValueError(f"File does not exist: {filename}")
        
        # 根据文件扩展名确定MS级别
        if filename.lower().endswith('.ms1'):
            level = 1
        elif filename.lower().endswith('.ms2'):
            level = 2
        else:
            raise ValueError(f"Unsupported file format: {filename}")
        
        ms_obj = MSFileObject(level=level)
        current_spectrum = None
        
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        for line in tqdm(lines, desc=f"Reading MS{level} file"):
            line = line.strip()
            
            # 跳过空行
            if not line:
                continue
            
            # 处理头信息行
            if line.startswith('H'):
                match = re.match(r'H\s+(\S+)\s+(.*)', line)
                if match:
                    key, value = match.groups()
                    ms_obj.set_metadata(key, value)
                continue
            
            # 开始新的谱图
            if line.startswith('S'):
                # 保存之前的谱图
                if current_spectrum:
                    ms_obj.add_spectrum(current_spectrum)
                
                # 解析S行
                parts = line.split()
                if len(parts) >= 3:
                    current_spectrum = MSSpectrum(level=level)
                    current_spectrum.scan_number = int(parts[2])
                    
                    # 如果是MS2，还需要解析前体离子m/z
                    if level == 2 and len(parts) >= 4:
                        current_spectrum.precursor_mz = float(parts[3])
                continue
            
            # 处理信息行
            if line.startswith('I') and current_spectrum:
                match = re.match(r'I\s+(\S+)\s+(.*)', line)
                if match:
                    key, value = match.groups()
                    if key == "RTime":
                        current_spectrum.retention_time = float(value)
                    else:
                        current_spectrum.set_additional_info(key, value)
                continue
            
            # 处理电荷行（仅MS2）
            if line.startswith('Z') and current_spectrum and level == 2:
                parts = line.split()
                if len(parts) >= 2:
                    current_spectrum.precursor_charge = int(parts[1])
                continue
            
            # 处理峰值数据
            if current_spectrum:
                try:
                    parts = line.split()
                    if len(parts) >= 2:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        current_spectrum.add_peak(mz, intensity)
                except ValueError:
                    pass
        
        # 添加最后一个谱图
        if current_spectrum:
            ms_obj.add_spectrum(current_spectrum)
        
        return ms_obj
    
    def read_to_msobjects(self, filename):
        """
        读取MS1/MS2文件并转换为MSObject对象列表
        
        Args:
            filename: MS1/MS2文件路径
            
        Returns:
            list: MSObject对象列表
        """
        # 读取MS文件
        ms_file_obj = self.read(filename)
        
        # 转换为MSObject列表
        ms_objects = []
        for spectrum in tqdm(ms_file_obj.spectra, desc="Converting to MSObjects"):
            ms_obj = MSFileConverter.to_msobject(spectrum)
            ms_objects.append(ms_obj)
        
        return ms_objects 