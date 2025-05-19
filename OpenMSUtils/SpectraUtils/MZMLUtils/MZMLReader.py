from lxml import etree
import os
import multiprocessing as mp
from tqdm import tqdm
from .MZMLObject import MZMLObject, Spectrum
import concurrent.futures
from .MZMLConverter import MZMLConverter

class MZMLReader(object):
    def __init__(self):
        super().__init__()

    def read(self, filename, parse_spectra=True, parallel=False, num_processes=None):
        """
        读取MZML文件并解析为MZMLObject对象
        
        Args:
            filename: mzML文件路径
            parse_spectra: 是否解析spectra列表，默认为True
            parallel: 是否使用并行处理解析spectra，默认为False
            num_processes: 并行处理的进程数，默认为None（使用CPU核心数）
            
        Returns:
            MZMLObject: 包含mzML数据的对象
        """
        xml_tree = etree.parse(filename)
        root = xml_tree.getroot()
        
        # 处理命名空间
        if root.tag.endswith('indexedmzML'):
            # 获取mzML节点
            for child in root:
                if child.tag.endswith('mzML'):
                    mzml_root = child
                    break
            # 创建MZMLObject，但不解析spectra
            if not parse_spectra:
                return MZMLObject(mzml_root, parse_spectra=False)
            
            # 如果需要解析spectra
            if parallel:
                # 使用并行处理解析spectra
                mzml_obj = MZMLObject(mzml_root, parse_spectra=False)
                self._parse_spectra_parallel(filename, mzml_obj, root, num_processes)
                return mzml_obj
            else:
                return MZMLObject(mzml_root, parse_spectra=True)
        else:
            # 直接创建MZMLObject并解析
            return MZMLObject(root, parse_spectra=parse_spectra)
    
    def read_to_msobjects(self, filename, parallel=False, num_processes=None):
        """
        读取MZML文件并解析为MSObject对象列表
        
        Args:
            filename: mzML文件路径
            parallel: 是否使用并行处理解析spectra，默认为False
            num_processes: 并行处理的进程数，默认为None（使用CPU核心数）
            
        Returns:
            list: MSObject对象列表
        """
        # 先读取为MZMLObject
        mzml_obj = self.read(filename, parse_spectra=True, parallel=parallel, num_processes=num_processes)
        
        # 将Spectrum对象转换为MSObject
        ms_objects = []
        if mzml_obj.run and mzml_obj.run.spectra_list:
            for spectrum in tqdm(mzml_obj.run.spectra_list, desc="Converting to MSObjects"):
                ms_obj = MZMLConverter.to_msobject(spectrum)
                ms_objects.append(ms_obj)
        
        return ms_objects

    def _parse_spectra_parallel(self, filename, mzml_obj, root, num_processes=None):
        """
        并行解析spectra
        
        Args:
            filename: mzML文件路径
            mzml_obj: MZMLObject对象
            root: XML根节点
            num_processes: 并行处理的进程数，默认为None（使用CPU核心数）
        """
        # 检查是否为indexedmzML
        if root.tag.endswith('indexedmzML'):
            # 使用索引并行解析spectra
            index_list, end_offset = self._get_offset_list(root)
            
            if num_processes is None:
                num_processes = mp.cpu_count()
            
            # 将索引分成多个块
            chunk_size = max(1, len(index_list) // num_processes)
            chunks = [index_list[i:i + chunk_size] for i in range(0, len(index_list), chunk_size)]
            
            # 创建线程池
            with concurrent.futures.ThreadPoolExecutor(max_workers=num_processes) as executor:
                # 使用map方法并行处理每个块
                results = list(tqdm(
                    executor.map(
                        lambda chunk: self._parse_spectra_chunk(filename, chunk, end_offset),
                        chunks
                    ),
                    total=len(chunks),
                    desc="Processing chunks"
                ))
                
                # 收集结果
                all_spectra = []
                for result in results:
                    all_spectra.extend(result)
                
                # 使用属性访问器设置spectra_list
                if mzml_obj.run:
                    mzml_obj.run.spectra_list = all_spectra
        else:
            # 如果不是indexedmzML，使用XML元素并行解析
            run_elem = None
            for child in root:
                if child.tag.endswith('run'):
                    run_elem = child
                    break
            if run_elem is not None:
                spectrum_list_elem = None
                for child in run_elem:
                    if child.tag.endswith('spectrumList'):
                        spectrum_list_elem = child
                        break
                if spectrum_list_elem is not None:
                    spectrum_elems = [child for child in spectrum_list_elem if child.tag.endswith('spectrum')]
                    
                    if num_processes is None:
                        num_processes = mp.cpu_count()
                    
                    # 将spectrum元素分成多个块
                    chunk_size = max(1, len(spectrum_elems) // num_processes)
                    chunks = [spectrum_elems[i:i + chunk_size] for i in range(0, len(spectrum_elems), chunk_size)]

                    # 创建线程池
                    with concurrent.futures.ThreadPoolExecutor(max_workers=num_processes) as executor:
                        # 并行处理每个块
                        results = list(tqdm(
                            executor.map(
                                lambda chunk: self._parse_spectrum_elems(chunk),
                                chunks
                            ),
                            total=len(chunks),
                            desc="Processing spectrum elements"
                        ))
                        
                        # 收集结果
                        all_spectra = []
                        for result in results:
                            all_spectra.extend(result)
                        
                        # 使用属性访问器设置spectra_list
                        if mzml_obj.run:
                            mzml_obj.run.spectra_list = all_spectra

    def _parse_spectrum_elems(self, spectrum_elems):
        """
        解析spectrum元素列表
        
        Args:
            spectrum_elems: spectrum元素列表
            
        Returns:
            list: Spectrum对象列表
        """
        spectra = []
        for spectrum_elem in spectrum_elems:
            try:
                spectrum = Spectrum(spectrum_elem)
                spectra.append(spectrum)
            except Exception as e:
                print(f"Error parsing spectrum element: {e}")
        
        return spectra

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
        
        index_list_elem = None
        for child in root:
            if child.tag.endswith('indexList'):
                index_list_elem = child
                break
        if index_list_elem is None:
            raise ValueError("No indexList found in the root element")
        
        # 获取spectrum索引
        spectrum_index_elem = None
        for child in index_list_elem:
            if child.tag.endswith('index') and child.get('name') == 'spectrum':
                spectrum_index_elem = child
                break
        if spectrum_index_elem is None:
            raise ValueError("No spectrum index found in the indexList element")
        
        # 遍历所有offset节点
        for offset_elem in [child for child in spectrum_index_elem if child.tag.endswith('offset')]:
            offset_value = int(offset_elem.text)
            id_ref = offset_elem.get('idRef')
            offset_list.append({
                'idRef': id_ref,
                'offset': offset_value
            })
        
        # 获取文件结束偏移量
        chromatogram_index_elem = None
        for child in index_list_elem:
            if child.tag.endswith('index') and child.get('name') == 'chromatogram':
                chromatogram_index_elem = child
                break
        if chromatogram_index_elem is not None:
            offset_elems = None
            for child in chromatogram_index_elem:
                if child.tag.endswith('offset'):
                    offset_elems = child
                    break
            if offset_elems is not None:
                end_offset = int(offset_elems.text) + 10000  # 添加一个足够大的值

        # 如果仍然没有找到结束偏移量，使用文件大小
        if end_offset is None:
            end_offset = os.path.getsize(root.base)
        
        return offset_list, end_offset

    def _parse_spectra_chunk(self, filename, offset_chunk, end_offset):
        """
        解析一个spectra块
        
        Args:
            filename: mzML文件路径
            offset_chunk: 偏移量块
            end_offset: 结束偏移量
            
        Returns:
            list: Spectrum对象列表
        """
        spectra = []
        with open(filename, 'rb') as file:
            for i, offset_info in enumerate(offset_chunk):
                file.seek(offset_info['offset'])
                
                # 确定读取的长度
                if i < len(offset_chunk) - 1:
                    read_length = offset_chunk[i+1]['offset'] - offset_info['offset']
                else:
                    read_length = end_offset - offset_info['offset']
                
                # 读取数据
                data = file.read(read_length)
                
                # 提取spectrum XML
                spectrum_start = data.find(b'<spectrum')
                spectrum_end = data.find(b'</spectrum>') + len(b'</spectrum>')
                if spectrum_start >= 0 and spectrum_end > 0:
                    spectrum_data = data[spectrum_start:spectrum_end]
                    
                    # 解析为Spectrum对象
                    try:
                        spectrum_elem = etree.fromstring(spectrum_data)
                        spectrum = Spectrum(spectrum_elem)
                        spectra.append(spectrum)
                    except Exception as e:
                        print(f"Error parsing spectrum at offset {offset_info['offset']}: {e}")
        
        return spectra

if __name__ == "__main__":
    file_path = 'D:\\code\\Python\\MS\\NADataFormer\\rawData\\20181121a_HAP1_tRNA_19.mzML'
    reader = MZMLReader()
    ms_data = reader.read(file_path)

    print(ms_data)
