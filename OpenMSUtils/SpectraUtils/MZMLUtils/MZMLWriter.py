from lxml import etree
import os
from tqdm import tqdm
from .MZMLObject import MZMLObject, Spectrum, Run
from .MZMLConverter import MZMLConverter

class MZMLWriter(object):
    def __init__(self):
        super().__init__()
    
    def write(self, mzml_obj, filename, write_index=True):
        """
        将MZMLObject写入mzML文件
        
        Args:
            mzml_obj: MZMLObject对象
            filename: 输出文件路径
            write_index: 是否写入索引，默认为True
            
        Returns:
            bool: 写入是否成功
        """
        try:
            # 创建mzML根元素
            mzml_root = etree.Element("mzML", nsmap=mzml_obj.nsmap)
            for key, value in mzml_obj.attrib.items():
                mzml_root.set(key, value)
            
            # 添加CV列表
            if mzml_obj.cv_list is not None:
                mzml_root.append(mzml_obj.cv_list)
            
            # 添加文件描述
            if mzml_obj.file_description is not None:
                mzml_root.append(mzml_obj.file_description)
            
            # 添加可引用参数组列表
            if mzml_obj.referenceable_param_group_list is not None:
                mzml_root.append(mzml_obj.referenceable_param_group_list)
            
            # 添加样本列表
            if mzml_obj.sample_list is not None:
                mzml_root.append(mzml_obj.sample_list)
            
            # 添加仪器配置列表
            if mzml_obj.instrument_configuration_list is not None:
                mzml_root.append(mzml_obj.instrument_configuration_list)
            
            # 添加软件列表
            if mzml_obj.software_list is not None:
                mzml_root.append(mzml_obj.software_list)
            
            # 添加数据处理列表
            if mzml_obj.data_processing_list is not None:
                mzml_root.append(mzml_obj.data_processing_list)
            
            # 添加采集列表
            if mzml_obj.acquisition_list is not None:
                mzml_root.append(mzml_obj.acquisition_list)
            
            # 添加运行信息
            if mzml_obj.run is not None:
                run_elem = mzml_obj.run.to_xml()
                mzml_root.append(run_elem)
            
            # 创建XML树
            tree = etree.ElementTree(mzml_root)
            
            if write_index:
                # 创建indexedmzML根元素
                indexed_root = etree.Element("indexedmzML", nsmap=mzml_obj.nsmap)
                indexed_root.append(mzml_root)
                
                # 写入文件（不包含索引）
                temp_filename = filename + ".temp"
                tree.write(temp_filename, pretty_print=True, xml_declaration=True, encoding="utf-8")
                
                # 读取文件并创建索引
                offsets = self._create_index(temp_filename)
                
                # 创建索引列表
                index_list = etree.SubElement(indexed_root, "indexList", count=str(len(offsets)))
                
                # 添加spectrum索引
                spectrum_index = etree.SubElement(index_list, "index", name="spectrum")
                for offset_info in offsets:
                    offset_elem = etree.SubElement(spectrum_index, "offset", idRef=offset_info['idRef'])
                    offset_elem.text = str(offset_info['offset'])
                
                # 添加文件校验和（可选）
                fileChecksum = etree.SubElement(indexed_root, "fileChecksum")
                fileChecksum.text = "0"
                
                # 写入最终文件
                indexed_tree = etree.ElementTree(indexed_root)
                indexed_tree.write(filename, pretty_print=True, xml_declaration=True, encoding="utf-8")
                
                # 删除临时文件
                os.remove(temp_filename)
            else:
                # 直接写入文件
                tree.write(filename, pretty_print=True, xml_declaration=True, encoding="utf-8")
            
            return True
        except Exception as e:
            print(f"Error writing mzML file: {e}")
            return False
    
    def write_from_msobjects(self, ms_objects, filename, metadata=None, write_index=True):
        """
        从MSObject列表创建并写入mzML文件
        
        Args:
            ms_objects: MSObject对象列表
            filename: 输出文件路径
            metadata: 元数据字典，包含文件描述、仪器配置等信息
            write_index: 是否写入索引，默认为True
            
        Returns:
            bool: 写入是否成功
        """
        # 创建基本的MZMLObject
        mzml_obj = MZMLObject()
        
        mzml_obj.nsmap = {
        None: "http://psi.hupo.org/ms/mzml",
        "xsi": "http://www.w3.org/2001/XMLSchema-instance",
        "schemaLocation": "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"
        }
        mzml_obj.attrib = {
            "version": "1.1.0",
            "xmlns": "http://psi.hupo.org/ms/mzml",
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xsi:schemaLocation": "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"
        }
        # 添加元数据（如果提供）
        if metadata:
            if 'file_description' in metadata:
                mzml_obj.file_description = metadata['file_description']
            if 'referenceable_param_group_list' in metadata:
                mzml_obj.referenceable_param_group_list = metadata['referenceable_param_group_list']
            if 'sample_list' in metadata:
                mzml_obj.sample_list = metadata['sample_list']
            if 'instrument_configuration_list' in metadata:
                mzml_obj.instrument_configuration_list = metadata['instrument_configuration_list']
            if 'software_list' in metadata:
                mzml_obj.software_list = metadata['software_list']
            if 'data_processing_list' in metadata:
                mzml_obj.data_processing_list = metadata['data_processing_list']
        
        # 创建Run对象
        run = Run()
        run.attrib = {"id": "run1", "defaultInstrumentConfigurationRef": "IC1"}
        
        # 将MSObject转换为Spectrum并添加到Run
        spectra = []
        for ms_obj in tqdm(ms_objects, desc="Converting MSObjects to Spectra"):
            spectrum = MZMLConverter.from_msobject(ms_obj)
            spectra.append(spectrum)
        
        run.spectra_list = spectra
        mzml_obj.run = run
        
        # 写入文件
        return self.write(mzml_obj, filename, write_index)
    
    def _create_index(self, filename):
        """
        创建mzML文件的索引
        
        Args:
            filename: mzML文件路径
            
        Returns:
            list: 包含偏移量信息的列表
        """
        offsets = []
        
        with open(filename, 'rb') as file:
            content = file.read()
            
            # 查找所有spectrum标签的位置
            start_pos = 0
            while True:
                spectrum_start = content.find(b'<spectrum ', start_pos)
                if spectrum_start == -1:
                    break
                
                # 提取spectrum ID
                id_start = content.find(b'id="', spectrum_start) + 4
                id_end = content.find(b'"', id_start)
                spectrum_id = content[id_start:id_end].decode('utf-8')
                
                offsets.append({
                    'idRef': spectrum_id,
                    'offset': spectrum_start
                })
                
                start_pos = spectrum_start + 10
        
        return offsets 