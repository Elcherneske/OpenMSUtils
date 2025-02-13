from lxml import etree
import re
import base64
import struct
import zlib

from .MSObject import MSObject

class MZMLReader(object):
    def __init__(self):
        super().__init__()

    def read(self, filename):
        """
        读取MZML文件并解析为MSObject对象列表
        """
        meta_data = {}
        spectrum = []
        xml_tree = etree.parse(filename)
        root = xml_tree.getroot()
        if root.tag.endswith('indexedmzML'): # 如果文件中存在indexList标签
            index_list, end_offset = self._get_offset_list(root)
            with open(filename, 'rb') as file:
                for i in range(len(index_list)):
                    file.seek(index_list[i]['offset'])
                    # 读取从当前offset到下一个offset之间的内容
                    if i == len(index_list) - 1:
                        data = file.read(end_offset - index_list[i]['offset']).split(b'</spectrum>')[0] + b'</spectrum>\n'
                    else:
                        data = file.read(index_list[i+1]['offset'] - index_list[i]['offset'])
                    ms_obj = self._parse_spectrum_to_msobject(data)
                    spectrum.append(ms_obj)
        
        return meta_data, spectrum

    def _get_offset_list(self, root):
        """
        从XML根节点获取所有offset值并构建列表
        Args:
            root: XML根节点
        Returns:
            list: 包含所有offset值的列表
        """
        offset_list = []
        for child in root:
            if child.tag.endswith('indexList'):
                index_root = child[0]
                break
        if index_root is None:
            raise ValueError("No indexList found in the root element")
        
        # 遍历所有offset节点
        for offset_elem in index_root:
            if offset_elem.tag.endswith('offset'):
                offset_value = offset_elem.text
                id_ref = offset_elem.get('idRef')
                offset_list.append({
                    'idRef': id_ref,
                    'offset': int(offset_value)
                })

        index_root = child[1]  # 获取chromatogram_index节点
        # 遍历所有offset节点, usually one
        for offset_elem in index_root:
            if offset_elem.tag.endswith('offset'):
                end_offset = offset_elem.text
                break

        if end_offset is None:
            raise ValueError("No end offset found in the indexList element")
                
        return offset_list, int(end_offset)

    def _parse_spectrum_to_msobject(self, spectrum_str):
        """
        将XML格式的质谱数据解析为MSObject对象
        Args:
            spectrum_str: XML格式的质谱字符串
        Returns:
            MSObject: 包含质谱数据的MSObject对象
        """
        # 解析XML字符串
        spectrum_elem = etree.fromstring(spectrum_str)
        
        # 创建MSObject对象
        ms_obj = MSObject()
        
        # 获取scan number
        id_str = spectrum_elem.get('id')
        scan_match = re.search(r'scan=(\d+)', id_str)
        if scan_match:
            scan_number = int(scan_match.group(1))
        else:
            # 如果没有匹配到scan=xxx格式,则为-1
            scan_number = -1
        ms_obj.set_scan_number(scan_number)
        
        # 遍历cvParam节点获取基本信息
        for cv_param in spectrum_elem.findall(".//cvParam"):
            name = cv_param.get('name')
            value = cv_param.get('value', '')
            
            # 获取MS level
            if name == "ms level":
                ms_obj.set_level(int(value))
                
            # 获取retention time
            elif name == "scan start time":
                ms_obj.set_retention_time(float(value))
                
            # 存储其他信息到additional_info
            else:
                ms_obj.add_additional_info(name, value)
        
        # 获取precursor信息
        precursorList = spectrum_elem.find(".//precursorList")
        if precursorList is not None:
            if len(precursorList) > 1:
                print("Warning: Multiple precursorList found in the spectrum element")
            for precursor in precursorList:
                selected_ion_list = precursor.find(".//selectedIonList")
                if selected_ion_list is not None:
                    for selected_ion in selected_ion_list:
                        mz = None
                        charge = None
                        intensity = None
                        for cv_param in selected_ion.findall("cvParam"):
                            if cv_param.get('name') == "selected ion m/z":
                                mz = float(cv_param.get('value'))
                            elif cv_param.get('name') == "charge state":
                                charge = int(cv_param.get('value'))
                            elif cv_param.get('name') == "peak intensity":
                                intensity = float(cv_param.get('value'))
                        if mz is not None:
                            if charge is not None:
                                ms_obj.set_precursor_ion(mz=mz, charge=charge)
                            else:
                                ms_obj.set_precursor_ion(mz=mz, charge=0)
                                
        # 获取mz和intensity数组
        binary_list = spectrum_elem.find("binaryDataArrayList")
        if binary_list is not None:
            binary_arrays = binary_list.findall("binaryDataArray")
            if len(binary_arrays) >= 2:
                # 解码base64编码的数据
                # 获取m/z和intensity数组的编码格式和压缩信息
                
                mz_cv_params = binary_arrays[0].findall("cvParam")
                mz_is_64_bit = any(cv_param.get('name') == "64-bit float" for cv_param in mz_cv_params)
                mz_is_compressed = any(cv_param.get('name') == "zlib compression" for cv_param in mz_cv_params)

                # 解码binary数据
                mz_binary_data = base64.b64decode(binary_arrays[0].find("binary").text)
                if mz_is_compressed:
                    mz_binary_data = zlib.decompress(mz_binary_data)

                # 根据位数解包数据
                num_floats = len(mz_binary_data) // (8 if mz_is_64_bit else 4)
                if mz_is_64_bit:
                    mz_values = struct.unpack('d' * num_floats, mz_binary_data)
                else:
                    mz_values = struct.unpack('f' * num_floats, mz_binary_data)

                # 获取intensity数组
                intensity_cv_params = binary_arrays[1].findall("cvParam")
                intensity_is_64_bit = any(cv_param.get('name') == "64-bit float" for cv_param in intensity_cv_params)
                intensity_is_compressed = any(cv_param.get('name') == "zlib compression" for cv_param in intensity_cv_params)
                
                # 解码binary数据
                intensity_binary_data = base64.b64decode(binary_arrays[1].find("binary").text)
                if intensity_is_compressed:
                    intensity_binary_data = zlib.decompress(intensity_binary_data)

                # 根据位数解包数据
                num_floats = len(intensity_binary_data) // (8 if intensity_is_64_bit else 4)
                if intensity_is_64_bit:
                    intensity_values = struct.unpack('d' * num_floats, intensity_binary_data)
                else:
                    intensity_values = struct.unpack('f' * num_floats, intensity_binary_data)

                # 添加峰值数据
                for mz, intensity in zip(mz_values, intensity_values):
                    ms_obj.add_peak(mz=mz, intensity=intensity)
        
        return ms_obj

if __name__ == "__main__":
    file_path = 'D:\\code\\Python\\MS\\NADataFormer\\rawData\\20181121a_HAP1_tRNA_19.mzML'
    reader = MZMLReader()
    ms_data = reader.read(file_path)

    print(ms_data)
