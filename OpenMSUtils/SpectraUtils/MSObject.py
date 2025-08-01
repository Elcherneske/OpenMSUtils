class MSObject:
    def __init__(self):
        self._level = 1
        self._scan = {
            'scan_number': -1,
            'rt': 0.0,
            'dt': 0.0,
            'scan_window': (0.0, 0.0),
        }
        self._precursor = None
        self._peaks = []
    
    @property
    def level(self):
        return self._level

    @property
    def peaks(self):
        return self._peaks

    @property
    def precursor(self):
        return self._precursor

    @property
    def scan(self):
        return self._scan

    @property
    def scan_number(self):
        if self._scan is not None and 'scan_number' in self._scan:
            return self._scan['scan_number']
        else:
            return None

    @property
    def retention_time(self):
        if self._scan is not None and 'rt' in self._scan:
            return self._scan['rt']
        else:
            return None
    
    @property
    def drift_time(self):
        if self._scan is not None and 'dt' in self._scan:
            return self._scan['dt']
        else:
            return None
    
    @property
    def precursor_mz(self):
        if self._precursor is not None and 'mz' in self._precursor:
            return self._precursor['mz']
        else:
            return None
    
    @property
    def precursor_charge(self):
        if self._precursor is not None and 'charge' in self._precursor:
            return self._precursor['charge']
        else:
            return None
    
    @property
    def precursor_window(self):
        if self._precursor is not None and 'isolation_window' in self._precursor:
            return self._precursor['isolation_window']
        else:
            return None

    def add_peak(self, mz:float, intensity:float):
        self._peaks.append((mz, intensity))
    
    def clear_peaks(self):
        self._peaks = []
    
    def set_peaks(self, peaks:list[tuple[float, float, int]] | list[tuple[float, float]]):
        self._peaks = peaks
    
    def sort_peaks(self):
        self._peaks.sort(key=lambda x: x[0])
    
    def set_level(self, level:int):
        self._level = level
    
    def set_precursor(
            self, 
            mz:float, 
            charge:int, 
            ref_scan_number:int, 
            activation_method:str, 
            activation_energy:float, 
            isolation_window:tuple[float, float]
        ):
        self._precursor = {
            'mz': mz,
            'charge': charge,
            'ref_scan_number': ref_scan_number,
            'activation_method': activation_method,
            'activation_energy': activation_energy,
            'isolation_window': isolation_window,
        }

    def set_scan(
            self, 
            scan_number:int, 
            retention_time:float, 
            drift_time:float, 
            scan_window:tuple[float, float]
        ):
        self._scan['scan_number'] = scan_number
        self._scan['rt'] = retention_time
        self._scan['dt'] = drift_time
        self._scan['scan_window'] = scan_window
    
if __name__ == "__main__":
    pass

