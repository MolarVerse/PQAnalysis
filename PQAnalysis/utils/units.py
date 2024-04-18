from unum import Unum
from unum.units import *

cal = Unum.unit('cal', J/4.184)
kcal = Unum.unit('kcal', cal*1e3)

mole = 6.02214076e23

cal_per_mole = Unum.unit('cal', J/4.184/mol)
kcal_per_mole = Unum.unit('kcal', cal_per_mole*1e3)
