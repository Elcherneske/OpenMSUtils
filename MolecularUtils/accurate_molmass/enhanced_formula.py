from functools import cached_property
from decimal import Decimal
from .molmass import Formula
from .elements import ELECTRON, ELEMENTS, Isotope

class EnhancedFormula(Formula):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    @cached_property
    def mass(self) -> float:
        """Average relative molecular mass.

        The sum of the relative atomic masses of all atoms and charges in
        the formula.
        Equals the molar mass in g/mol, that is, the mass of one mole of
        substance.

        >>> Formula('H').mass
        1.007941
        >>> Formula('H+').mass
        1.007392...
        >>> Formula('SO4_2-').mass
        96.06351...
        >>> Formula('12C').mass
        12.0
        >>> Formula('C8H10N4O2').mass
        194.1909...
        >>> Formula('C48H32AgCuO12P2Ru4').mass
        1438.404...

        """
        result = Decimal("0.0")
        for symbol, massnumber_counts in self._elements.items():
            ele = ELEMENTS[symbol]
            for massnumber, count in massnumber_counts.items():
                if massnumber:
                    result += Decimal(f"{ele.isotopes[massnumber].mass * count}")
                else:
                    result += Decimal(f"{ele.mass * count}")
        return float(result - Decimal(f"{ELECTRON.mass * self._charge}"))
    
    @cached_property
    def isotope(self) -> Isotope:
        """Isotope composed of most abundant elemental isotopes.

        >>> Formula('C').isotope.mass
        12.0
        >>> Formula('13C').isotope.massnumber
        13
        >>> Formula('C48H32AgCuO12P2Ru4').isotope
        Isotope(mass=1439.588..., abundance=0.00205..., massnumber=1440...)

        """
        result = Isotope(Decimal(f"{-ELECTRON.mass * self._charge}"), Decimal("1.0"), 0, self._charge)
        for symbol, massnumber_counts in self._elements.items():
            ele = ELEMENTS[symbol]
            for massnumber, count in massnumber_counts.items():
                if massnumber != 0:
                    isotope = ele.isotopes[massnumber]
                else:
                    isotope = ele.isotopes[ele.nominalmass]
                result.mass += Decimal(f"{isotope.mass * count}")
                result.massnumber += isotope.massnumber * count
                result.abundance *= Decimal(f"{isotope.abundance ** count}")
        result.mass = float(result.mass)
        result.abundance = float(result.abundance)
        return result