from tqdm import tqdm
import os
import multiprocessing
from .MSFileConverter import MSFileConverter
from .MGFConverter import MGFConverter
from .MZMLConverter import MZMLConverter

def mzml_chunk_to_spectra_objects(chunk):
    return [MZMLConverter.to_spectra_object(spectrum) for spectrum in chunk]

def mgf_chunk_to_spectra_objects(chunk):
    return [MGFConverter.to_spectra_object(spectrum) for spectrum in chunk]

def msfile_chunk_to_spectra_objects(chunk):
    return [MSFileConverter.to_spectra_object(spectrum) for spectrum in chunk]

class MSReader(object):
    def __init__(self, thread_num=None):
        super().__init__()
        if thread_num is None:
            self.thread_num = multiprocessing.cpu_count()
        else:   
            self.thread_num = thread_num

    def read_to_spectra_objects(self, filename):
        """
        读取质谱文件并解析为SpectraObject对象
        
        Args:
            filename: 质谱文件路径
            
        Returns:
            SpectraObject: 包含Spectra数据的对象列表
        """
        if not os.path.exists(filename):
            raise ValueError(f"File does not exist: {filename}")
        
        if filename.endswith('.mzML'):
            return self._mzml_to_spectra_objects(filename)
        elif filename.endswith('.mgf'):
            return self._mgf_to_spectra_objects(filename)
        elif filename.endswith('.ms1') or filename.endswith('.ms2'):
            return self._msfile_to_spectra_objects(filename)
        else:
            raise ValueError(f"Unsupported file type: {filename}")
    
    def _mzml_to_spectra_objects(self, filename):
        """
        读取MZML文件并解析为SpectraObject对象列表

        Args:
            filename: mzML文件路径

        Returns:
            list: SpectraObject对象列表
        """
        from lxml import etree

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
                results = list(tqdm(pool.imap(mzml_chunk_to_spectra_objects, spectrum_chunks), total=len(spectrum_chunks), desc="Converting to MSObjects"))
                for ms_list in results:
                    ms_objects.extend(ms_list)
        else:
            ms_objects = [MZMLConverter.to_spectra_object(spectrum) for spectrum in spectrum_list]

        return ms_objects
    
    def _mgf_to_spectra_objects(self, filename):
        """
        读取MGF文件并解析为SpectraObject对象
        
        Args:
            filename: MGF文件路径
            
        Returns:
            SpectraObject: 包含Spectra数据的对象
        """
        if not os.path.exists(filename):
            raise ValueError(f"File does not exist: {filename}")
        
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        spectra_indices = [i for i, line in enumerate(lines) if line.startswith('BEGIN IONS')]
        spectra_indices.append(len(lines))

        spectras = []
        for i in range(len(spectra_indices) - 1):
            start_index = spectra_indices[i]
            end_index = spectra_indices[i + 1]
            spectras.append(lines[start_index:end_index])
        
        ms_objects = []
        if self.thread_num > 1:
            chunk_size = (len(spectras) + self.thread_num - 1) // self.thread_num
            spectras_chunks = [spectras[i*chunk_size:(i+1)*chunk_size] for i in range(self.thread_num)]

            with multiprocessing.Pool(processes=self.thread_num) as pool:
                results = list(tqdm(pool.imap(mgf_chunk_to_spectra_objects, spectras_chunks), total=len(spectras_chunks), desc="Converting lines to MSObject"))
                for ms_list in results:
                    ms_objects.extend(ms_list)
        else:
            ms_objects = [MGFConverter.to_spectra_object(spectrum) for spectrum in spectras]

        return ms_objects
    
    def _msfile_to_spectra_objects(self, filename):
        """
        读取MS1/MS2文件并解析为SpectraObject对象
        
        Args:
            filename: MS1/MS2文件路径
            
        Returns:
            SpectraObject: 包含Spectra数据的对象
        """
        if not os.path.exists(filename):
            raise ValueError(f"File does not exist: {filename}")
        
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        spectra_indices = [i for i, line in enumerate(lines) if line.startswith('S')]
        spectra_indices.append(len(lines))

        spectras = []
        for i in range(len(spectra_indices) - 1):
            start_index = spectra_indices[i]
            end_index = spectra_indices[i + 1]
            spectras.append(lines[start_index:end_index])
        
        ms_objects = []
        if self.thread_num > 1:
            chunk_size = (len(spectras) + self.thread_num - 1) // self.thread_num
            spectras_chunks = [spectras[i*chunk_size:(i+1)*chunk_size] for i in range(self.thread_num)]

            with multiprocessing.Pool(processes=self.thread_num) as pool:
                results = list(tqdm(pool.imap(msfile_chunk_to_spectra_objects, spectras_chunks), total=len(spectras_chunks), desc="Converting lines to MSObject"))
                for ms_list in results:
                    ms_objects.extend(ms_list)
        else:
            ms_objects = [MSFileConverter.to_spectra_object(spectrum) for spectrum in spectras]

        return ms_objects
