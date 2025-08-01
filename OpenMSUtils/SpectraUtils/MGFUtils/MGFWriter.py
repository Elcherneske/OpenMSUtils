from .MGFConverter import MGFConverter

class MGFWriter(object):
    def __init__(self):
        super().__init__()
    
    def write_from_msobjects(self, ms_objects, filename):
        """
        从MSObject列表创建并写入MGF文件
        
        Args:
            ms_objects: MSObject对象列表
            filename: 输出文件路径
            
        Returns:
            bool: 写入是否成功
        """
        lines = []
        for ms_obj in ms_objects:
            lines.extend(MGFConverter.from_msobject(ms_obj))
        
        with open(filename, 'w') as file:
            file.write('\n'.join(lines))
        
        return True