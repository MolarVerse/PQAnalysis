"""
A module containing different units related to the unum subpackage.
"""

from unum import Unum
from unum.units import J

cal = Unum.unit('cal', J/4.184)
kcal = Unum.unit('kcal', cal*1e3)

mole = 6.02214076e23

cal_per_mole = Unum.unit('cal/mol', J/4.184/mole)
kcal_per_mole = Unum.unit('kcal/mol', cal_per_mole*1e3)
