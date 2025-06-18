# OpenMSUtils

This open-source repository is used for easily reading and processing MS files and Fasta files.

## Installation

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

## Interfaces

### SpectraUtils

This package provides various classes for reading and processing MS files, including MZML, MGF, MS1, and MS2.

1. MZMLReader: read MZML file

```python
read (file_path, parse_spectra=True, parallel=False, num_processes=None) -> MZMLObject
read_to_msobjects (file_path, parallel=False, num_processes=None) -> list[MSObject]
```

2. MGFReader: read MGF file

```python
read (file_path) -> MGFObject
read_to_msobjects (file_path) -> list[MSObject]
```

3. MSFileReader: read MS1/MS2 file

```python
read (file_path) -> MSFileObject
read_to_msobjects (file_path) -> list[MSObject]
```

4. SpectraPlotter: plot spectra

```python
plot_mz (ms_object: MSObject)
plot_ion_mobility (ion_mobility_spectrum: dict, time_range: tuple[float, float] = None, time_bins: int = 500, mz_range: tuple[float, float] = None, mz_bins: int = 500)
plot_xics (precursor_xics: List[XICResult], fragment_xics: List[XICResult], output_file: Optional[str] = None)
```

5. XICSExtractor: extract XICs

```python
extract_xics (df: pd.DataFrame) -> List[tuple[List[XICResult], List[XICResult]]]
```

### FastaUtils

This package provides various classes for reading and processing Fasta files.

1. FastaReader: read Fasta file

```python
read (file_path) -> FastaObject
```

2. FastaWriter: write Fasta file

```python
write (sequences: dict, file_path: str)
```

### AnalysisUtils

This package provides various classes for analyzing MS data.

1. FDRUtils: calculate FDR

```python
calculate_fdr (score, label, target_fdr=0.01, top_n=20) -> int, float
```

### MolecularUtils

This package provides various classes for processing molecular data.

1. NucleicAcidUtils: process nucleic acid sequences

```python
class Oligonucleotide:
__init__ (sequence: str, deoxidation=False)
set_end_modifications (end_3_modification: tuple[str, str], end_5_modification: tuple[str, str])
add_modification (index: int, modification: Modification)
set_charge (charge: int)
set_adduct (adduct: str)
set_fragments_type (fragments_type: list[str])
mass -> float
mz -> float
fragments -> list[tuple[str, float]]
```

2. ProteinUtils: process protein sequences

```python
class Peptide:
__init__ (sequence: str)
set_end_modifications (end_C_modification: tuple[str, str], end_N_modification: tuple[str, str])
add_modification (index: int, modification: Modification)
set_charge (charge: int)
set_adduct (adduct: str)
set_fragments_type (fragments_type: list[str])
mass -> float
mz -> float
fragments -> list[tuple[str, float]]
```

3. ModificationUtils: process modifications

```python
class Modification:
__init__ (name: str, formula: str)
mass -> float

class ModificationUtils:
parse_modification_file (file_path: str) -> pd.DataFrame
find_modifications_by_mass (modifications: pd.DataFrame, mass: float, tolerance: float = 0.0001) -> List[Modification]
parse_modified_sequence (modified_sequence: str) -> Tuple[str, Dict[int, str]]
format_modified_sequence (sequence: str, modifications: Dict[int, str]) -> str
```

4. DecoyUtils: generate decoy sequences

```python
generate_decoy (sequence: str, modifications: Optional[Dict[int, str]] = None) -> Tuple[str, Dict[int, str]]
generate_decoy_batch (sequences: List[str], modifications_list: Optional[List[Dict[int, str]]] = None) -> List[Tuple[str, Dict[int, str]]]
calculate_similarity (seq1: str, seq2: str) -> float
```





