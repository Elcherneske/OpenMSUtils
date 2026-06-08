import base64
import re
import struct
import zlib
from lxml import etree
from .MSObject import SpectraObject

class MZMLConverter:
    @staticmethod
    def to_spectra_object(spectrum: etree._Element | str) -> SpectraObject:
        """
        将mzML的Spectrum对象转换为SpectraObject
        
        Args:
            spectrum: MZMLObject中的Spectrum对象
            
        Returns:
            SpectraObject对象
        """
        # 创建MSObject
        spectra_object = SpectraObject()

        if not isinstance(spectrum, etree._Element):
            spectrum = etree.fromstring(spectrum)

        # 获取 namespace
        nsmap = spectrum.nsmap
        if None in nsmap:
            ns = {'ns': nsmap[None]}
        else:
            # 取第一个命名空间
            prefix = next(iter(nsmap))
            ns = {prefix: nsmap[prefix]}

        # 解析 scan number
        scan_number = -1
        if 'index' in spectrum.attrib:
            scan_number = int(spectrum.attrib.get('index')) + 1
        if 'id' in spectrum.attrib:
            id_str = spectrum.attrib.get('id', '')
            if 'scan=' in id_str:
                scan_number_str = id_str.split('scan=')[1]
                match = re.search(r'(\d+)', scan_number_str)
                if match:
                    scan_number = int(match.group(1))   

        # 拆分子Element, 提取所有cvParam
        mz_array_element = None
        intensity_array_element = None
        cv_params = spectrum.findall('.//ns:cvParam', namespaces=ns)

        for child in spectrum.findall('.//ns:binaryDataArray', namespaces=ns):
            for cv_param in child.findall('.//ns:cvParam', namespaces=ns):
                if cv_param.attrib.get('accession') == 'MS:1000514':  # m/z array
                    mz_array_element = child
                    break
                if cv_param.attrib.get('accession') == 'MS:1000515':  # intensity array
                    intensity_array_element = child
                    break
        
        ms_level = 1
        retention_time = 0.0
        drift_time = 0.0
        scan_window = (0.0, 0.0)
        isolation_target = 0.0
        precursor_mz = 0.0
        precursor_charge = 0
        activation_method = 'unknown'
        activation_energy = 0.0
        isolation_window = (0.0, 0.0)

        for cv_param in cv_params:
            if cv_param.attrib.get('accession') == 'MS:1000511':  # ms level
                ms_level = int(cv_param.attrib.get('value', '0'))
            elif cv_param.attrib.get('accession') == 'MS:1000016':  # scan start time
                value = float(cv_param.attrib.get('value', '0'))
                unit = cv_param.attrib.get('unitAccession', '')
                # 转换为秒
                if unit == 'UO:0000031':  # minutes
                    retention_time = value * 60
                elif unit == 'UO:0000028':  # milliseconds
                    retention_time = value / 1000
                else:
                    retention_time = value
            elif cv_param.attrib.get('accession') == 'MS:1002476':  # ion mobility drift time
                value = float(cv_param.attrib.get('value', '0'))
                unit = cv_param.attrib.get('unitAccession', '')
                # 转换为秒
                if unit == 'UO:0000031':  # minutes
                    drift_time = value * 60
                elif unit == 'UO:0000028':  # milliseconds
                    drift_time = value / 1000
                else:
                    drift_time = value
            elif cv_param.attrib.get('accession') == 'MS:1000501':  # scan window lower limit
                scan_window = (float(cv_param.attrib.get('value', '0')), scan_window[1])
            elif cv_param.attrib.get('accession') == 'MS:1000500':  # scan window upper limit
                scan_window = (scan_window[0], float(cv_param.attrib.get('value', '0')))
            elif cv_param.attrib.get('accession') == 'MS:1000827':  # isolation window target m/z
                isolation_target = float(cv_param.attrib.get('value', '0'))
            elif cv_param.attrib.get('accession') == 'MS:1000828':  # isolation window lower offset
                isolation_window = (float(cv_param.attrib.get('value', '0')), isolation_window[1])
            elif cv_param.attrib.get('accession') == 'MS:1000829':  # isolation window upper offset
                isolation_window = (isolation_window[0], float(cv_param.attrib.get('value', '0')))
            elif cv_param.attrib.get('accession') == 'MS:1000744':  # selected ion m/z
                precursor_mz = float(cv_param.attrib.get('value', '0'))
            elif cv_param.attrib.get('accession') == 'MS:1000041':  # charge state
                precursor_charge = int(cv_param.attrib.get('value', '0'))
            elif cv_param.attrib.get('accession') in ['MS:1000133', 'MS:1000134', 'MS:1000422', 'MS:1000250']:
                activation_method = cv_param.attrib.get('name', 'unknown')
            elif cv_param.attrib.get('accession') == 'MS:1000045':  # collision energy
                activation_energy = float(cv_param.attrib.get('value', '0'))

        isolation_window = (isolation_target - isolation_window[0], isolation_target + isolation_window[1])

        # 设置MS级别
        spectra_object.set_level(ms_level)
        
        spectra_object.set_scan(scan_number, retention_time, drift_time, scan_window)

        ref_scan_number = -1
        precursor_element_list = spectrum.findall('.//ns:precursor', namespaces=ns)
        if len(precursor_element_list) > 0:
            precursor_element = precursor_element_list[0]
            if 'spectrumRef' in precursor_element.attrib:
                ref_id = precursor_element.attrib.get('spectrumRef', '')
                if 'scan=' in ref_id:
                    ref_scan_number = int(ref_id.split('scan=')[1].split()[0]) + 1

            spectra_object.set_precursor(precursor_mz, precursor_charge, ref_scan_number, activation_method, activation_energy, isolation_window)
            
        
        if mz_array_element is None:
            raise ValueError("m/z array is not found")
        
        if intensity_array_element is None:
            raise ValueError("intensity array is not found")

        mz_precision = 64  # 默认为双精度
        mz_compression = False
        for cv_param in mz_array_element.findall('.//ns:cvParam', namespaces=ns):
            if cv_param.attrib.get('accession') == 'MS:1000521':  # 32-bit float
                mz_precision = 32
            elif cv_param.attrib.get('accession') == 'MS:1000523':  # 64-bit float
                mz_precision = 64
            elif cv_param.attrib.get('accession') == 'MS:1000574':  # zlib compression
                mz_compression = True
        
        mz_array = mz_array_element.findall('.//ns:binary', namespaces=ns)[0].text

        intensity_precision = 64
        intensity_compression = False
        for cv_param in intensity_array_element.findall('.//ns:cvParam', namespaces=ns):
            if cv_param.attrib.get('accession') == 'MS:1000521':  # 32-bit float
                intensity_precision = 32
            elif cv_param.attrib.get('accession') == 'MS:1000523':  # 64-bit float
                intensity_precision = 64
            elif cv_param.attrib.get('accession') == 'MS:1000574':  # zlib compression
                intensity_compression = True
        
        intensity_array = intensity_array_element.findall('.//ns:binary', namespaces=ns)[0].text

        if mz_array is None or intensity_array is None:
            spectra_object.clear_peaks()
            return spectra_object

        mz_array = base64.b64decode(mz_array)
        if mz_compression:
            mz_array = zlib.decompress(mz_array)
        
        if mz_precision == 32:
            fmt = 'f' * (len(mz_array) // 4)
            mz_array = struct.unpack(fmt, mz_array)
        else:  # 64-bit
            fmt = 'd' * (len(mz_array) // 8)
            mz_array = struct.unpack(fmt, mz_array)

        intensity_array = base64.b64decode(intensity_array)
        if intensity_compression:
            intensity_array = zlib.decompress(intensity_array)
        
        if intensity_precision == 32:
            fmt = 'f' * (len(intensity_array) // 4)
            intensity_array = struct.unpack(fmt, intensity_array)
        else:  # 64-bit
            fmt = 'd' * (len(intensity_array) // 8)
            intensity_array = struct.unpack(fmt, intensity_array)
        
        # 如果找到了m/z和intensity数组，则添加峰值
        if mz_array and intensity_array:
            assert len(mz_array) == len(intensity_array)
            spectra_object.clear_peaks()  # 清除现有峰值
            peaks = [(mz, intensity) for mz, intensity in zip(mz_array, intensity_array)]
            peaks.sort(key=lambda x: x[0])
            spectra_object.set_peaks(peaks)
        
        return spectra_object
    
    @staticmethod
    def from_spectra_object(spectra_object: SpectraObject) -> etree._Element:
        """
        将SpectraObject转换为mzML的Spectrum对象
        
        Args:
            spectra_object: SpectraObject对象
            
        Returns:
            MZMLObject中的Spectrum对象
        """
        spectrum = etree.Element('spectrum')
        
        # 设置基本属性
        spectrum.attrib = {
            'index': str(spectra_object.scan_number - 1),
            'id': f'scan={spectra_object.scan_number}',
            'defaultArrayLength': str(len(spectra_object.peaks))
        }
        
        # 添加MS级别
        ms_level_param = etree.Element('cvParam')
        ms_level_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000511',
            'name': 'ms level',
            'value': str(spectra_object.level)
        }
        spectrum.append(ms_level_param)
        
        # 添加scan_list
        scan_list = etree.Element('scanList')
        scan_list.attrib = {
            'count': str(1)
        }

        scan = etree.Element('scan')
        rt_param = etree.Element('cvParam')
        rt_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000016',
            'name': 'scan start time',
            'value': str(spectra_object.retention_time / 60),  # 转换为分钟
            'unitCvRef': 'UO',
            'unitAccession': 'UO:0000031',
            'unitName': 'minute'
        }
        scan.append(rt_param)

        if spectra_object.drift_time > 0:
            dt_param = etree.Element('cvParam')
            dt_param.attrib = {
                'cvRef': 'MS',
                'accession': 'MS:1002476',
                'name': 'ion mobility drift time',
                'value': str(spectra_object.scan.drift_time * 1000),  # 转换为毫秒
                'unitCvRef': 'UO',
                'unitAccession': 'UO:0000028',
                'unitName': 'millisecond'
            }
            scan.append(dt_param)
        
        scan_window_list = etree.Element('scanWindowList')
        scan_window_list.attrib = {
            'count': str(1)
        }
        scan_window = etree.Element('scanWindow')
        low_param = etree.Element('cvParam')
        low_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000501',
            'name': 'scan window lower limit',
            'value': str(spectra_object.scan.scan_window[0])
        }
        high_param = etree.Element('cvParam')
        high_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000500',
            'name': 'scan window upper limit',
            'value': str(spectra_object.scan.scan_window[1])
        }
        scan_window.append(low_param)
        scan_window.append(high_param)
        scan_window_list.append(scan_window)
        scan.append(scan_window_list)
        scan_list.append(scan)

        # 添加precursor
        if not spectra_object.precursor is None:
            precursor_list = etree.Element('precursorList')
            precursor_list.attrib = {
                'count': str(1)
            }
            precursor = etree.Element('precursor')
            precursor.attrib = {
                'spectrumRef': f'scan={spectra_object.precursor.ref_scan_number}'
            }

            isolation_window = etree.Element('isolationWindow')
            target = (spectra_object.precursor.isolation_window[0] + spectra_object.precursor.isolation_window[1]) / 2
            low_offset = target - spectra_object.precursor.isolation_window[0]
            high_offset = spectra_object.precursor.isolation_window[1] - target
            target_param = etree.Element('cvParam')
            target_param.attrib = {
                'cvRef': 'MS',
                'accession': 'MS:1000827',
                'name': 'isolation window target m/z',
                'value': str(target)
            }
            isolation_window.append(target_param)
            low_param = etree.Element('cvParam')
            low_param.attrib = {
                'cvRef': 'MS',
                'accession': 'MS:1000828',
                'name': 'isolation window lower offset',
                'value': str(low_offset)
            }
            isolation_window.append(low_param)
            high_param = etree.Element('cvParam')
            high_param.attrib = {
                'cvRef': 'MS',
                'accession': 'MS:1000829',
                'name': 'isolation window upper offset',
                'value': str(high_offset)
            }
            isolation_window.append(high_param)
            precursor.append(isolation_window)
            
            selected_ion_list = etree.Element('selectedIonList')
            selected_ion_list.attrib = {
                'count': str(1)
            }
            selected_ion = etree.Element('selectedIon')
            mz_param = etree.Element('cvParam')
            mz_param.attrib = {
                'cvRef': 'MS',
                'accession': 'MS:1000744',
                'name': 'selected ion m/z',
                'value': str(spectra_object.precursor_mz)
            }
            selected_ion.append(mz_param)
            if spectra_object.precursor.charge != 0:
                charge_param = etree.Element('cvParam')
                charge_param.attrib = {
                    'cvRef': 'MS',
                    'accession': 'MS:1000041',
                    'name': 'charge state',
                    'value': str(spectra_object.precursor_charge)
                }
                selected_ion.append(charge_param)
            selected_ion_list.append(selected_ion)
            precursor.append(selected_ion_list)
            
            activation = etree.Element('activation')
            method_param = etree.Element('cvParam')
            method_param.attrib = {
                'cvRef': 'MS',
                'accession': 'MS:1000133',
                'name': 'collision energy',
                'value': str(spectra_object.precursor['activation_method'])
            }
            activation.append(method_param)
            if spectra_object.precursor['activation_energy'] > 0:
                energy_param = etree.Element('cvParam')
                energy_param.attrib = {
                    'cvRef': 'MS',
                    'accession': 'MS:1000045',
                    'name': 'collision energy',
                    'value': str(spectra_object.precursor['activation_energy'])
                }
                activation.append(energy_param)
            precursor.append(activation)
            precursor_list.append(precursor)
            spectrum.append(precursor_list)

        # 添加peak_list
        mz_values = [peak[0] for peak in spectra_object.peaks]
        intensity_values = [peak[1] for peak in spectra_object.peaks]
        binary_data_array_list = etree.Element('binaryDataArrayList')
        binary_data_array_list.attrib = {
            'count': str(2)
        }
        mz_array = etree.Element('binaryDataArray')
        array_type_param = etree.Element('cvParam')
        array_type_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000514',
            'name': 'm/z array'
        }
        mz_array.append(array_type_param)
        mz_precision_param = etree.Element('cvParam')
        mz_precision_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000523',
            'name': '64-bit float'
        }
        mz_array.append(mz_precision_param)
        mz_compression_param = etree.Element('cvParam')
        mz_compression_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000574',
            'name': 'zlib compression'
        }
        mz_array.append(mz_compression_param)

        mz_binary = struct.pack('d' * len(mz_values), *mz_values)
        mz_compressed = zlib.compress(mz_binary)
        mz_encoded = base64.b64encode(mz_compressed).decode('ascii')
        mz_array.append(etree.Element('binary', text=mz_encoded))
        binary_data_array_list.append(mz_array)

        intensity_array = etree.Element('binaryDataArray')
        array_type_param = etree.Element('cvParam')
        array_type_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000515',
            'name': 'intensity array'
        }
        intensity_array.append(array_type_param)
        intensity_precision_param = etree.Element('cvParam')
        intensity_precision_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000523',
            'name': '64-bit float'
        }
        intensity_array.append(intensity_precision_param)
        intensity_compression_param = etree.Element('cvParam')
        intensity_compression_param.attrib = {
            'cvRef': 'MS',
            'accession': 'MS:1000574',
            'name': 'zlib compression'
        }
        intensity_array.append(intensity_compression_param)

        intensity_binary = struct.pack('d' * len(intensity_values), *intensity_values)
        intensity_compressed = zlib.compress(intensity_binary)
        intensity_encoded = base64.b64encode(intensity_compressed).decode('ascii')
        intensity_array.append(etree.Element('binary', text=intensity_encoded))
        binary_data_array_list.append(intensity_array)

        spectrum.append(binary_data_array_list)
        return spectrum