import os
import multiprocessing
from tqdm import tqdm
from .MGFConverter import MGFConverter

def process_chunk(chunk):
    return [MGFConverter.to_msobject(spectrum) for spectrum in chunk]

class MGFReader(object):
    def __init__(self, thread_num=None):
        super().__init__()
        if thread_num is None:
            self.thread_num = multiprocessing.cpu_count()
        else:
            self.thread_num = thread_num

    def read_to_msobjects(self, filename):
        """
        读取MS1/MS2文件并解析为MSFileObject对象
        
        Args:
            filename: MS1/MS2文件路径
            
        Returns:
            MSFileObject: 包含MS数据的对象
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
                results = list(tqdm(pool.imap(process_chunk, spectras_chunks), total=len(spectras_chunks), desc="Converting lines to MSObject"))
                for ms_list in results:
                    ms_objects.extend(ms_list)
        else:
            ms_objects = [MGFConverter.to_msobject(spectrum) for spectrum in spectras]

        return ms_objects