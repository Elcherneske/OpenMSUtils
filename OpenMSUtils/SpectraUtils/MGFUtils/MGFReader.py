import os
from tqdm import tqdm
from .MGFObject import MGFObject, MGFSpectrum
from .MGFConverter import MGFConverter

class MGFReader(object):
    def __init__(self):
        super().__init__()
    
    def read(self, filename):
        """
        读取MGF文件并解析为MGFObject对象
        
        Args:
            filename: MGF文件路径
            
        Returns:
            MGFObject: 包含MGF数据的对象
        """
        if not os.path.exists(filename) or not filename.lower().endswith('.mgf'):
            raise ValueError(f"Invalid file name: {filename}")
        
        mgf_obj = MGFObject()
        current_spectrum = None
        
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        for line in tqdm(lines, desc="Reading MGF file"):
            line = line.strip()
            
            # 跳过空行
            if not line:
                continue
            
            # 处理注释行（元数据）
            if line.startswith('#'):
                if '=' in line:
                    key, value = line[1:].strip().split('=', 1)
                    mgf_obj.set_metadata(key.strip(), value.strip())
                continue
            
            # 开始新的谱图
            if line == "BEGIN IONS":
                current_spectrum = MGFSpectrum()
                continue
            
            # 结束当前谱图
            if line == "END IONS":
                if current_spectrum:
                    mgf_obj.add_spectrum(current_spectrum)
                    current_spectrum = None
                continue
            
            # 处理谱图参数
            if '=' in line and current_spectrum:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                
                if key == "TITLE":
                    current_spectrum.title = value
                elif key == "PEPMASS":
                    # 处理可能包含强度的PEPMASS
                    pepmass_parts = value.split()
                    current_spectrum.pepmass = float(pepmass_parts[0])
                elif key == "CHARGE":
                    # 处理电荷格式，如"2+"或"3-"
                    if value.endswith('+'):
                        current_spectrum.charge = int(value[:-1])
                    elif value.endswith('-'):
                        current_spectrum.charge = -int(value[:-1])
                    else:
                        try:
                            current_spectrum.charge = int(value)
                        except ValueError:
                            pass
                elif key == "RTINSECONDS":
                    current_spectrum.rtinseconds = float(value)
                else:
                    current_spectrum.set_additional_info(key, value)
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
        
        return mgf_obj
    
    def read_to_msobjects(self, filename):
        """
        读取MGF文件并转换为MSObject对象列表
        
        Args:
            filename: MGF文件路径
            
        Returns:
            list: MSObject对象列表
        """
        # 读取MGF文件
        mgf_obj = self.read(filename)
        
        # 转换为MSObject列表
        ms_objects = []
        for spectrum in tqdm(mgf_obj.spectra, desc="Converting to MSObjects"):
            ms_obj = MGFConverter.to_msobject(spectrum)
            ms_objects.append(ms_obj)
        
        return ms_objects 