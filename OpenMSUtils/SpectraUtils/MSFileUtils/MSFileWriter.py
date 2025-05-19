import os
from tqdm import tqdm
from .MSFileObject import MSFileObject, MSSpectrum
from .MSFileConverter import MSFileConverter

class MSFileWriter(object):
    def __init__(self):
        super().__init__()
    
    def write(self, ms_file_obj, filename):
        """
        将MSFileObject写入MS1/MS2文件
        
        Args:
            ms_file_obj: MSFileObject对象
            filename: 输出文件路径
            
        Returns:
            bool: 写入是否成功
        """
        try:
            # 检查文件扩展名是否与MS级别匹配
            if ms_file_obj.level == 1 and not filename.lower().endswith('.ms1'):
                filename = filename + '.ms1'
            elif ms_file_obj.level == 2 and not filename.lower().endswith('.ms2'):
                filename = filename + '.ms2'
            
            # 转换为MS格式字符串
            ms_string = ms_file_obj.to_ms_string()
            
            # 写入文件
            with open(filename, 'w') as file:
                file.write(ms_string)
            
            return True
        except Exception as e:
            print(f"Error writing MS file: {e}")
            return False
    
    def write_from_msobjects(self, ms_objects, filename, metadata=None):
        """
        从MSObject列表创建并写入MS1/MS2文件
        
        Args:
            ms_objects: MSObject对象列表
            filename: 输出文件路径
            metadata: 元数据字典，可选
            
        Returns:
            bool: 写入是否成功
        """
        if not ms_objects:
            raise ValueError("No MS objects provided")
        
        # 确定MS级别
        level = ms_objects[0].level
        
        # 检查所有对象是否具有相同的级别
        for ms_obj in ms_objects:
            if ms_obj.level != level:
                raise ValueError("All MS objects must have the same level")
        
        # 检查文件扩展名是否与MS级别匹配
        if level == 1 and not filename.lower().endswith('.ms1'):
            filename = filename + '.ms1'
        elif level == 2 and not filename.lower().endswith('.ms2'):
            filename = filename + '.ms2'
        
        # 创建MSFileObject
        ms_file_obj = MSFileObject(level=level)
        
        # 添加元数据
        if metadata:
            for key, value in metadata.items():
                ms_file_obj.set_metadata(key, value)
        
        # 将MSObject转换为MSSpectrum并添加到MSFileObject
        for ms_obj in tqdm(ms_objects, desc=f"Converting MSObjects to MS{level}"):
            ms_spectrum = MSFileConverter.from_msobject(ms_obj)
            ms_file_obj.add_spectrum(ms_spectrum)
        
        # 写入文件
        return self.write(ms_file_obj, filename) 