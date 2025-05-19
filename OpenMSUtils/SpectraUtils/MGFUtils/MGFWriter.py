from tqdm import tqdm
from .MGFObject import MGFObject, MGFSpectrum
from .MGFConverter import MGFConverter

class MGFWriter(object):
    def __init__(self):
        super().__init__()
    
    def write(self, mgf_obj, filename):
        """
        将MGFObject写入MGF文件
        
        Args:
            mgf_obj: MGFObject对象
            filename: 输出文件路径
            
        Returns:
            bool: 写入是否成功
        """
        try:
            # 转换为MGF格式字符串
            mgf_string = mgf_obj.to_mgf_string()
            
            # 写入文件
            with open(filename, 'w') as file:
                file.write(mgf_string)
            
            return True
        except Exception as e:
            print(f"Error writing MGF file: {e}")
            return False
    
    def write_from_msobjects(self, ms_objects, filename, metadata=None):
        """
        从MSObject列表创建并写入MGF文件
        
        Args:
            ms_objects: MSObject对象列表
            filename: 输出文件路径
            metadata: 元数据字典，可选
            
        Returns:
            bool: 写入是否成功
        """
        # 创建MGFObject
        mgf_obj = MGFObject()
        
        # 添加元数据
        if metadata:
            for key, value in metadata.items():
                mgf_obj.set_metadata(key, value)
        
        # 将MSObject转换为MGFSpectrum并添加到MGFObject
        for ms_obj in tqdm(ms_objects, desc="Converting MSObjects to MGF"):
            mgf_spectrum = MGFConverter.from_msobject(ms_obj)
            mgf_obj.add_spectrum(mgf_spectrum)
        
        # 写入文件
        return self.write(mgf_obj, filename) 