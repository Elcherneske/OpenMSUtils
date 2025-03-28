from .accurate_molmass import EnhancedFormula
from .ModificationUtils import Modification

class Nucleotide():
    def __init__(self, character, deoxidation=False):
        self.character = character
        self.deoxidation = deoxidation

    @property
    def mass(self):
        formula = EnhancedFormula("C5H8O6P")
        if self.character == "A":
            formula += EnhancedFormula("C5H4N5")
        elif self.character == "G":
            formula += EnhancedFormula("C5H4N5O")
        elif self.character == "C":
            formula += EnhancedFormula("C4H4N3O")
        elif self.character == "T":
            formula += EnhancedFormula("C5H5N2O2")
        elif self.character == "U":
            formula += EnhancedFormula("C4H5N2O2")

        if self.deoxidation:
            formula -= EnhancedFormula("O")

        return formula.isotope.mass
    
    @property
    def nucleobase(self):
        if self.character == "A":
            return EnhancedFormula("C5H4N5").isotope
        elif self.character == "G":
            return EnhancedFormula("C5H4N5O").isotope
        elif self.character == "C":
            return EnhancedFormula("C4H4N3O").isotope
        elif self.character == "T":
            return EnhancedFormula("C5H5N2O2").isotope
        elif self.character == "U":
            return EnhancedFormula("C4H5N2O2").isotope

class Oligonucleotide():
    """
    Oligonucleotide class
    """
    def __init__(self, sequence: str, deoxidation=False):
        """
        sequence: str (5' -> 3')
        """
        self._sequence = sequence
        self._deoxidation = deoxidation
        self._modifications = {}
        self._charge = 0
        self._adduct_mass = None
        self._end_3_modification = None
        self._end_5_modification = None
        self._fragments_type = None
    
    def set_end_modifications(self, end_3_modification: tuple[str, str], end_5_modification: tuple[str, str]):
        end_3_name, end_3_formula = end_3_modification
        end_5_name, end_5_formula = end_5_modification

        self._end_3_modification = Modification(end_3_name, end_3_formula)
        self._end_5_modification = Modification(end_5_name, end_5_formula)
    
    def add_modification(self, index: int, modification: Modification):
        self._modifications[index] = modification
    
    def set_charge(self, charge: int):
        self._charge = charge
    
    def set_adduct(self, adduct: str):
        if adduct == "@H+":
            self._adduct_mass = -EnhancedFormula("H+").isotope.mass
        elif adduct == "H+":
            self._adduct_mass = EnhancedFormula("H+").isotope.mass
        elif adduct == "Na+":
            self._adduct_mass = EnhancedFormula("Na+").isotope.mass
        elif adduct == "K+":
            self._adduct_mass = EnhancedFormula("K+").isotope.mass
        elif adduct == "NH4+":
            self._adduct_mass = EnhancedFormula("NH4+").isotope.mass
        else:
            raise ValueError(f"Invalid adduct: {adduct}")
    
    def set_fragments_type(self, fragments_type: list[str]):
        for fragment_type in fragments_type:
            if fragment_type not in ["a-B", "a", "b", "c", "d", "w", "x", "y", "z"]:
                raise ValueError(f"Invalid fragment type: {fragment_type}")
        self._fragments_type = fragments_type
    
    @property
    def mass(self):
        mass = self._end_5_modification.mass if self._end_5_modification else 0.0
        for i, nucleotide in enumerate(self._sequence):
            mass += Nucleotide(nucleotide, self._deoxidation).mass
            if i in self._modifications:
                mass += self._modifications[i].mass
        mass += self._end_3_modification.mass if self._end_3_modification else 0.0
        return mass

    @property
    def mz(self):
        if self._charge == 0:
            return self.mass
        else:
            if self._adduct_mass is None and self._charge < 0:
                self._adduct_mass = -EnhancedFormula("H+").isotope.mass
            elif self._adduct_mass is None and self._charge > 0:
                self._adduct_mass = EnhancedFormula("H+").isotope.mass
            return (self.mass + abs(self._charge) * self._adduct_mass)/ abs(self._charge)
    
    @property
    def fragments(self):
        adduct_mass = self._adduct_mass if self._adduct_mass is not None else -EnhancedFormula("H+").isotope.mass
        charge = self._charge if self._charge is not None else -1
        if not (adduct_mass < 0 and charge < 0) or (adduct_mass > 0 and charge > 0):
            raise ValueError(f"Invalid adduct/charge combination: adduct={adduct_mass}, charge={charge}")
        fragments_type = self._fragments_type if self._fragments_type is not None else ["a", "b", "c", "d", "w", "x", "y", "z"]
        charge_sign = 1 if charge > 0 else -1
        max_fragment_charge = abs(charge)

        fragments = {}
        for ion_type in fragments_type:
            fragment_masses = self._generate_fragments(ion_type)
            if max_fragment_charge == 1:
                max_fragment_charge = 2
            for z in range(1, max_fragment_charge):
                charge_symbol = '+' if charge_sign > 0 else '-'
                fragment_key = f"{ion_type}{z}{charge_symbol}"
                # Calculate m/z values for each fragment mass
                fragments[fragment_key] = [(mass + z * adduct_mass) / z for mass in fragment_masses]

        return fragments
    
    def _generate_fragments(self, ion_type: str):
        ion_mod = {
            'a': -EnhancedFormula("HPO3").isotope.mass, 
            'a_B': -EnhancedFormula("HPO3").isotope.mass, 
            'b': -EnhancedFormula("HPO2").isotope.mass,
            'c': 0.0,
            'd': EnhancedFormula("O").isotope.mass,
            'w': EnhancedFormula("HPO3").isotope.mass, 
            'x': EnhancedFormula("HPO2").isotope.mass, 
            'y': 0.0, 
            'z': -EnhancedFormula("O").isotope.mass
        }

        fragments = []
        if ion_type in ["a", "a_B", "b", "c", "d"]:
            fragment_mass = self._end_5_modification.mass if self._end_5_modification else 0.0
            for i, nucleotide in enumerate(self._sequence):
                fragment_mass += Nucleotide(nucleotide, self._deoxidation).mass
                if i in self._modifications:
                    fragment_mass += self._modifications[i].mass
                # Skip record the first and last nucleotide
                if i == 0 or i == len(self._sequence) - 1:
                    continue
                fragments.append(fragment_mass + ion_mod[ion_type])

        elif ion_type in ["w", "x", "y", "z"]:
            fragment_mass = self._end_3_modification.mass if self._end_3_modification else 0.0
            for i in reversed(range(len(self._sequence))):
                nucleotide = self._sequence[i]
                fragment_mass += Nucleotide(nucleotide, self._deoxidation).mass
                if i in self._modifications:
                    fragment_mass += self._modifications[i].mass
                # Skip record the first and last nucleotide
                if i == 0 or i == len(self._sequence) - 1:
                    continue
                fragments.append(fragment_mass + ion_mod[ion_type])

        return fragments
                

if __name__ == "__main__":
    pass

