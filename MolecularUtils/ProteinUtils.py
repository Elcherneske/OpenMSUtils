from .accurate_molmass import EnhancedFormula
from .ModificationUtils import Modification

amino_acids = {
    "A": {"formula": "C3H5NO", "mass": EnhancedFormula("C3H5NO").isotope.mass}, #Alanine
    "C": {"formula": "C3H5NOS", "mass": EnhancedFormula("C3H5NOS").isotope.mass}, #Cysteine
    "D": {"formula": "C4H5NO3", "mass": EnhancedFormula("C4H5NO3").isotope.mass}, #Aspartic acid
    "E": {"formula": "C5H7NO3", "mass": EnhancedFormula("C5H7NO3").isotope.mass}, #Glutamic acid
    "F": {"formula": "C9H9NO", "mass": EnhancedFormula("C9H9NO").isotope.mass}, #Phenylalanine
    "G": {"formula": "C2H3NO", "mass": EnhancedFormula("C2H3NO").isotope.mass}, #Glycine
    "H": {"formula": "C6H7N3O", "mass": EnhancedFormula("C6H7N3O").isotope.mass}, #Histidine
    "I": {"formula": "C6H11NO", "mass": EnhancedFormula("C6H11NO").isotope.mass}, #Isoleucine
    "K": {"formula": "C6H12N2O", "mass": EnhancedFormula("C6H12N2O").isotope.mass}, #Lysine
    "L": {"formula": "C6H11NO", "mass": EnhancedFormula("C6H11NO").isotope.mass}, #Leucine
    "M": {"formula": "C5H9NOS", "mass": EnhancedFormula("C5H9NOS").isotope.mass}, #Methionine
    "N": {"formula": "C4H6N2O2", "mass": EnhancedFormula("C4H6N2O2").isotope.mass}, #Asparagine
    "P": {"formula": "C5H7NO", "mass": EnhancedFormula("C5H7NO").isotope.mass}, #Proline
    "Q": {"formula": "C5H8N2O2", "mass": EnhancedFormula("C5H8N2O2").isotope.mass}, #Glutamine
    "R": {"formula": "C6H12N4O", "mass": EnhancedFormula("C6H12N4O").isotope.mass}, #Arginine
    "S": {"formula": "C3H5NO2", "mass": EnhancedFormula("C3H5NO2").isotope.mass}, #Serine
    "T": {"formula": "C4H7NO2", "mass": EnhancedFormula("C4H7NO2").isotope.mass}, #Threonine
    "U": {"formula": "C3H5NOSe", "mass": EnhancedFormula("C3H5NOSe").isotope.mass}, #selenocysteine
    "V": {"formula": "C5H9NO", "mass": EnhancedFormula("C5H9NO").isotope.mass}, #Valine
    "W": {"formula": "C11H10N2O", "mass": EnhancedFormula("C11H10N2O").isotope.mass}, #Tryptophan
    "Y": {"formula": "C9H9NO2", "mass": EnhancedFormula("C9H9NO2").isotope.mass}, #Tyrosine
}

ion_mod = {
    'a': {"formula": "[M-CO]", "mass": -EnhancedFormula("CO").isotope.mass}, 
    'b': {"formula": "[M]", "mass": 0.0},
    'c': {"formula": "[M+NH]", "mass": EnhancedFormula("NH").isotope.mass},
    'x': {"formula": "[M+CO]", "mass": EnhancedFormula("CO").isotope.mass}, 
    'y': {"formula": "[M]", "mass": 0.0}, 
    'z': {"formula": "[M-NH]", "mass": -EnhancedFormula("NH").isotope.mass}
}

class Peptide():
    def __init__(self, 
        sequence: str, 
        modifications: dict = None,
        charge: int = 0,
        adduct: str = '[M+H+]',
        fragments_type: list[str] = None,
        end_C_modification: Modification = None,
        end_N_modification: Modification = None
    ):
        """
        sequence: str (N-terminus -> C-terminus)
        """
        self._sequence = sequence
        self._charge = charge
        self._adduct = Modification('adduct', adduct)
        if modifications is None:
           self._modifications = {}
        else:
            self._modifications = modifications

        if end_C_modification is None:
            self._end_C_modification = Modification('C-term', '[M+OH]')
        else:
            self._end_C_modification = end_C_modification

        if end_N_modification is None:
            self._end_N_modification = Modification('N-term', '[M+H]')
        else:
            self._end_N_modification = end_N_modification

        if fragments_type is None:
            self._fragments_type = ['b', 'y']
        else:
            self._fragments_type = fragments_type
    
    def set_C_modification(self, end_C_modification: Modification):
        self._end_C_modification = end_C_modification
    
    def set_N_modification(self, end_N_modification: Modification):
        self._end_N_modification = end_N_modification
    
    def add_modification(self, index: int, modification: Modification):
        if index in self._modifications:
            raise ValueError(f"Modification at index {index} already exists")
        self._modifications[index] = modification
    
    def set_charge(self, charge: int):
        self._charge = charge
    
    def set_adduct(self, adduct: str):
        self._adduct = Modification('adduct', adduct)
    
    def set_fragments_type(self, fragments_type: list[str]):
        for fragment_type in fragments_type:
            if fragment_type not in ["a", "b", "c", "x", "y", "z"]:
                raise ValueError(f"Invalid fragment type: {fragment_type}")
        self._fragments_type = fragments_type

    @property
    def mass(self):
        total_mass = sum(amino_acids[aa]["mass"] for aa in self._sequence)
        for index, modification in self._modifications.items():
            total_mass += modification.mass
        total_mass += self._end_C_modification.mass
        total_mass += self._end_N_modification.mass
        return total_mass

    @property
    def mz(self):
        if self._charge == 0:
            raise ValueError("Charge is 0, no m/z can be calculated")
        return (self.mass + (abs(self._charge) / abs(self._adduct.charge)) * self._adduct.mass) / abs(self._charge)

    @property
    def fragments(self):
        if self._charge == 0:
            raise ValueError("Charge is 0, no fragments can be calculated")
            
        charge_sign = 1 if self._charge > 0 else -1
        max_fragment_charge = abs(self._charge)
        if max_fragment_charge == 1:
            max_fragment_charge = 2

        fragments = {}
        for ion_type in self._fragments_type:
            fragment_masses = self._generate_fragments(ion_type)
            for z in range(1, max_fragment_charge):
                charge_symbol = '+' if charge_sign > 0 else '-'
                fragment_key = f"{ion_type}{z}{charge_symbol}"
                fragments[fragment_key] = [(mass + (z / abs(self._adduct.charge)) * self._adduct.mass) / z for mass in fragment_masses]

        return fragments
    
    def _generate_fragments(self, ion_type: str):
        fragments = []
        if ion_type in ["a", "b", "c"]:
            fragment_mass = self._end_N_modification.mass
            for i, aa in enumerate(self._sequence):
                fragment_mass += amino_acids[aa]["mass"]
                if i in self._modifications:
                    fragment_mass += self._modifications[i].mass
                fragments.append(fragment_mass + ion_mod[ion_type]["mass"])

        if ion_type in ["x", "y", "z"]:
            fragment_mass = self._end_C_modification.mass
            for i in reversed(range(len(self._sequence))):
                aa = self._sequence[i]
                fragment_mass += amino_acids[aa]["mass"]
                if i in self._modifications:
                    fragment_mass += self._modifications[i].mass
                fragments.append(fragment_mass + ion_mod[ion_type]["mass"])

        return fragments

if __name__ == "__main__":
    pass
