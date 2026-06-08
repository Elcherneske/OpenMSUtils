This open-source repository is used for easily reading and processing MS files, proteins and nucleic acids.

# Installation

use pip to install the package

```bash
pip install git+https://github.com/Elcherneske/OpenMSUtils.git
```

or

```bash
git clone https://github.com/Elcherneske/OpenMSUtils.git
cd OpenMSUtils
pip install .
```

# Interfaces

## MSUtils
This module provides class for reading and writing MS files.
1. SpectraObject: a class for representing single spectrum data.

Main properties of this class:

- level: (read-only) The MS level of the spectrum (e.g., MS1, MS2).
- scan: (read-only) A dictionary containing scan information, including scan_number (scan id), rt (retention time), dt (drift time), and scan_window (tuple).
- precursor: (read-only) A dictionary containing precursor ion information, including mz (mass/charge), charge, ref_scan_number, activation_method, activation_energy, and isolation_window.
- peaks: (read-only) A numpy array containing the peak list, usually as (mz, intensity).
- scan_number: (read-only) The scan number, equivalent to scan['scan_number'].
- retention_time: (read-only) The retention time, equivalent to scan['rt'].
- drift_time: (read-only) The drift time, equivalent to scan['dt'].
- precursor_mz: (read-only) The m/z value of the precursor ion.
- precursor_charge: (read-only) The charge of the precursor ion.
- precursor_window: (read-only) The isolation window of the precursor ion.

2. MSWriter: a class for writing MS files.

init parameters:
- thread_num: the number of threads to use for writing the file.

functions:
- write_from_spectra_objects(spectra_objects: list[SpectraObject], filename: str): write the spectra objects to a MS file, support mzML, mgf, ms1, ms2 files.

3. MSReader: a class for reading MS files.

init parameters:
- thread_num: the number of threads to use for reading the file.

functions:
- read_to_spectra_objects(filename: str): read the file and return a list of SpectraObject, support mzML, mgf, ms1, ms2 files.

## SpectraUtils

This package provides various classes for processing MS files.

1. XICSExtractor: extract XICs

init parameters:
- ppm_tolerance: (optional) the ppm tolerance for the XIC extraction, default is 25.0.
- rt_bin_size: (optional) the retention time bin size for the XIC extraction, default is 1.0.
- num_threads: (optional) the number of threads to use for the XIC extraction, default is 1.
- min_scans: (optional) the minimum number of scans to be included in the XIC, default is 5.
- peak_boundary: (optional) the peak boundary for the XIC extraction, default is 0.2.
- mode: (optional) the mode for the XIC extraction, 'rt_range' or 'scan_window', default is 'rt_range'.

functions:
- extract_xics(mzml_file: str, df: pd.DataFrame) -> List[tuple[List[XICResult], List[XICResult]]]: extract the XICs from the SpectraObject list.

## AnalysisUtils

This package provides various classes for analyzing various data.

1. FDRUtils: calculate FDR

functions:
- calculate_fdr(score, label, target_fdr=0.01, top_n=20) -> int, float: return the number of targets below threshold and the minimum score threshold.

2. FastaUtils: process Fasta files

functions:
- read(filename: str) -> dict: {header: sequence, ...}
- write(sequences: dict, filename: str): write the sequences to a Fasta file.

3. DecoyUtils: generate decoy sequences

functions:
- generate_decoy(sequence: str, modifications: Optional[Dict[int, str]] = None, method: str = "reverse", keep_terminals: bool = True, similarity_threshold: float = 1.0, max_attempts: int = 10) -> Tuple[str, Dict[int, str]]: generate a decoy sequence, return the decoy sequence and the modifications dictionary.

## MolecularUtils

This package provides various classes for simulating molecular data.

1. Modification: a class for representing a modification.

init parameters:
- name: (required) the name of the modification.
- formula: (required) the formula of the modification, the format should be '[M+x]' or '[M-x]', where x is the chemical formula of the modification like 'H2O' or 'NH3'.

properties:
- mass: (read-only) the mass of the modification.
- charge: (read-only) the charge of the modification.
- name: (read-only) the name of the modification.
- formula: (read-only) the formula of the modification.

2. ModificationUtils: a class for processing modifications.

functions:
- parse_modification_file(file_path: str) -> pd.DataFrame: parse the modification file and return a DataFrame. # in development
- find_modifications_by_mass(modifications: pd.DataFrame, mass: float, tolerance: float = 0.0001) -> List[Modification]: find the modifications by mass. # in development
- parse_modified_sequence(modified_sequence: str) -> Tuple[str, Dict[int, str]]: parse the modified sequence and return the sequence and the modifications, for example, 'PEPTIDE(UniMod:1)' will be parsed to ('PEPTIDE', {6: 'UniMod:1'}).
- format_modified_sequence(sequence: str, modifications: Dict[int, str]) -> str: format the modified sequence and return the sequence, for example, ('PEPTIDE', {6: 'UniMod:1'}) will be formatted to 'PEPTIDE(UniMod:1)'.

3. Peptide: a class for representing a peptide.

init parameters:
- sequence: (required) the sequence of the peptide.
- modifications: (optional) the modifications of the peptide, default is None.
- charge: (optional) the charge of the peptide, default is 0.
- adduct: (optional) the adduct of the peptide, default is '[M+H+]'.
- fragments_type: (optional) the fragments type of the peptide, default is ['b', 'y'].

properties:
- mass: (read-only) the mass of the peptide.
- mz: (read-only) the m/z value of the peptide.
- fragments: (read-only) the fragments of the peptide.

functions:
- add_modification(index: int, modification: Modification): add a modification to the peptide.
- set_charge(charge: int): set the charge of the peptide.
- set_adduct(adduct: str): set the adduct of the peptide.
- set_fragments_type(fragments_type: list[str]): set the fragments type of the peptide.

4. NucleicAcidUtils: a class for representing a nucleic acid. # in development




