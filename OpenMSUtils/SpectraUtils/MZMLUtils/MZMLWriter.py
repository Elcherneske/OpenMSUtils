from lxml import etree
import os
from tqdm import tqdm
from .MZMLConverter import MZMLConverter

class MZMLWriter(object):
    def __init__(self):
        super().__init__()

    def write_from_msobjects(self, ms_objects, filename):
        """
        从MSObject列表创建并写入mzML文件，并可选添加index

        Args:
            ms_objects: MSObject对象列表
            filename: 输出文件路径
            write_index: 是否写入索引（indexedmzML）

        Returns:
            bool: 写入是否成功
        """
        try:
            # 创建mzML根节点
            mzml_obj = etree.Element("mzML")
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

            # 创建Run对象
            run = etree.Element("run")
            spectrum_list = etree.Element("spectrumList")
            spectrum_list.attrib = {"count": str(len(ms_objects))}
            for ms_obj in tqdm(ms_objects, desc="Converting MSObjects to Spectra"):
                spectrum = MZMLConverter.from_msobject(ms_obj)
                spectrum_list.append(spectrum)
            run.append(spectrum_list)
            mzml_obj.append(run)

            # 创建indexedmzML根元素
            indexed_root = etree.Element("indexedmzML", nsmap=mzml_obj.nsmap)
            indexed_root.append(mzml_obj)

            # 写入临时mzML文件
            temp_filename = filename + ".temp"
            mzml_obj.write(temp_filename, pretty_print=True, xml_declaration=True, encoding="utf-8")

            # 读取文件并创建索引
            offsets = self._create_index(temp_filename)

            # 创建索引列表
            index_list = etree.Element("indexList", count=str(len(offsets)))
            spectrum_index = etree.Element("index", name="spectrum")
            for offset_info in offsets:
                offset_elem = etree.Element("offset")
                offset_elem.text = str(offset_info)
                spectrum_index.append(offset_elem)
            index_list.append(spectrum_index)
            indexed_root.append(index_list)

            # 写入最终文件
            indexed_tree = etree.ElementTree(indexed_root)
            indexed_tree.write(filename, pretty_print=True, xml_declaration=True, encoding="utf-8")

            # 删除临时文件
            import os
            os.remove(temp_filename)

            return True
        except Exception as e:
            print(f"写入mzML文件出错: {e}")
            return False

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

                offsets.append(spectrum_start)

                start_pos = spectrum_start + 10

        return offsets