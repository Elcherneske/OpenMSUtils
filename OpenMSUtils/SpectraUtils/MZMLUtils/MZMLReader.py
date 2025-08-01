from lxml import etree
import os
import multiprocessing
from tqdm import tqdm
import concurrent.futures
from .MZMLConverter import MZMLConverter

def process_chunk(chunk):
    return [MZMLConverter.to_msobject(spectrum) for spectrum in chunk]

class MZMLReader(object):
    def __init__(self, thread_num=None):
        super().__init__()
        if thread_num is None:
            self.thread_num = multiprocessing.cpu_count()
        else:
            self.thread_num = thread_num
    
    def read_to_msobjects(self, filename):
        """
        读取MZML文件并解析为MSObject对象列表

        Args:
            filename: mzML文件路径

        Returns:
            list: MSObject对象列表
        """
        root = etree.parse(filename).getroot()

        # 获取 namespace
        nsmap = root.nsmap
        if None in nsmap:
            ns = {'ns': nsmap[None]}
        else:
            prefix = next(iter(nsmap))
            ns = {prefix: nsmap[prefix]}

        spectrum_list = root.findall('.//ns:spectrum', namespaces=ns)

        if self.thread_num > 1:
            manager = multiprocessing.Manager()
            spectrum_str_list = manager.list([etree.tostring(spectrum, encoding='utf-8', pretty_print=False) for spectrum in spectrum_list])

            if not spectrum_list:
                raise ValueError("No spectrum found in the mzML file")

            # 将 spectrum_str_list 分为 self.thread_num 份
            chunk_size = (len(spectrum_str_list) + self.thread_num - 1) // self.thread_num
            spectrum_chunks = [spectrum_str_list[i*chunk_size:(i+1)*chunk_size] for i in range(self.thread_num)]

            ms_objects = []

            with multiprocessing.Pool(processes=self.thread_num) as pool:
                results = list(tqdm(pool.imap(process_chunk, spectrum_chunks), total=len(spectrum_chunks), desc="Converting to MSObjects"))
                for ms_list in results:
                    ms_objects.extend(ms_list)
        else:
            ms_objects = [MZMLConverter.to_msobject(spectrum) for spectrum in spectrum_list]

        return ms_objects

    def _get_offset_list(self, root):
        """
        从XML根节点获取所有offset值并构建列表
        Args:
            root: XML根节点
        Returns:
            list: 包含所有offset值的列表 
            int: 结束偏移量
        """
        offset_list = []
        end_offset = None

        nsmap = root.nsmap
        if None in nsmap:
            ns = {'ns': nsmap[None]}
        else:
            prefix = next(iter(nsmap))
            ns = {prefix: nsmap[prefix]}
        
        index_elem_list = root.findall('.//ns:index', namespaces=ns)
        
        # 获取spectrum索引
        spectrum_index_elem = None
        for index_elem in index_elem_list:
            if index_elem.get('name') == 'spectrum':
                spectrum_index_elem = index_elem
                break

        if spectrum_index_elem is None:
            raise ValueError("No spectrum index found in the indexList element")
        
        # 遍历所有offset节点
        for offset_elem in spectrum_index_elem:
            if not offset_elem.tag.endswith('offset'):
                continue

            offset_list.append(int(offset_elem.text))
        
        end_offset = offset_list[-1] + 10000  # 添加一个足够大的值

        return offset_list, end_offset

if __name__ == "__main__":
    pass
