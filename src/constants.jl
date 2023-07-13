# General constants

using Unitful, UnitfulAtomic

# definitions from RASPA

# const ANGSTROM = 1e-10
# const NANO_SECOND = 1e-9
# const PICO_SECOND = 1e-12
# const ATOMIC_MASS_UNIT = 1.6605402e-27
# const MOLAR_GAS_CONSTANT = 8.314464919
# const AVOGADRO_CONSTANT = 6.0221419947e23
# const ENERGY_CONVERSION_FACTOR = ATOMIC_MASS_UNIT * ANGSTROM^2/PICO_SECOND^2
# const ENERGY_TO_KELVIN = ENERGY_CONVERSION_FACTOR * AVOGADRO_CONSTANT / MOLAR_GAS_CONSTANT
# const ELECTRONIC_CHARGE_UNIT = 1.60217733e-19
# const ELECTRIC_CONSTANT = 8.85418781762039e-12
# const DielectricConstantOfTheMedium = 1.0
# const COULOMBIC_CONVERSION_FACTOR = ELECTRONIC_CHARGE_UNIT^2/(4π*ELECTRIC_CONSTANT*ANGSTROM*ENERGY_CONVERSION_FACTOR*DielectricConstantOfTheMedium)

const ENERGY_CONVERSION_FACTOR = 1u"u * Å^2 / ps^2"
const ENERGY_TO_KELVIN = NoUnits(ENERGY_CONVERSION_FACTOR/u"k_au*K")*u"K"
const VACUUM_PERMITTIVITY = 8.85418781762039e-12u"F/m"
const COULOMBIC_CONVERSION_FACTOR = NoUnits(1u"e_au^2"/(4π*VACUUM_PERMITTIVITY*u"Å"*ENERGY_CONVERSION_FACTOR))

# used for grid interpolation
const COEFF = Float64.(
   [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -3 3 0 0 0 0 0 0 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    2 -2 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    9 -9 -9 9 0 0 0 0 6 3 -6 -3 0 0 0 0 6 -6 3 -3 0 0 0 0 0 0 0 0 0 0 0 0 4 2 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -6 6 6 -6 0 0 0 0 -3 -3 3 3 0 0 0 0 -4 4 -2 2 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -6 6 6 -6 0 0 0 0 -4 -2 4 2 0 0 0 0 -3 3 -3 3 0 0 0 0 0 0 0 0 0 0 0 0 -2 -1 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    4 -4 -4 4 0 0 0 0 2 2 -2 -2 0 0 0 0 2 -2 2 -2 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 3 0 0 0 0 0 0 -2 -1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 -9 -9 9 0 0 0 0 0 0 0 0 0 0 0 0 6 3 -6 -3 0 0 0 0 6 -6 3 -3 0 0 0 0 4 2 2 1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6 6 6 -6 0 0 0 0 0 0 0 0 0 0 0 0 -3 -3 3 3 0 0 0 0 -4 4 -2 2 0 0 0 0 -2 -2 -1 -1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6 6 6 -6 0 0 0 0 0 0 0 0 0 0 0 0 -4 -2 4 2 0 0 0 0 -3 3 -3 3 0 0 0 0 -2 -1 -2 -1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 -4 -4 4 0 0 0 0 0 0 0 0 0 0 0 0 2 2 -2 -2 0 0 0 0 2 -2 2 -2 0 0 0 0 1 1 1 1 0 0 0 0;
    -3 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 -3 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    9 -9 0 0 -9 9 0 0 6 3 0 0 -6 -3 0 0 0 0 0 0 0 0 0 0 6 -6 0 0 3 -3 0 0 0 0 0 0 0 0 0 0 4 2 0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -6 6 0 0 6 -6 0 0 -3 -3 0 0 3 3 0 0 0 0 0 0 0 0 0 0 -4 4 0 0 -2 2 0 0 0 0 0 0 0 0 0 0 -2 -2 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 0 -1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 -9 0 0 -9 9 0 0 0 0 0 0 0 0 0 0 6 3 0 0 -6 -3 0 0 0 0 0 0 0 0 0 0 6 -6 0 0 3 -3 0 0 4 2 0 0 2 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6 6 0 0 6 -6 0 0 0 0 0 0 0 0 0 0 -3 -3 0 0 3 3 0 0 0 0 0 0 0 0 0 0 -4 4 0 0 -2 2 0 0 -2 -2 0 0 -1 -1 0 0;
    9 0 -9 0 -9 0 9 0 0 0 0 0 0 0 0 0 6 0 3 0 -6 0 -3 0 6 0 -6 0 3 0 -3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 2 0 2 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 9 0 -9 0 -9 0 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 0 3 0 -6 0 -3 0 6 0 -6 0 3 0 -3 0 0 0 0 0 0 0 0 0 4 0 2 0 2 0 1 0;
    -27 27 27 -27 27 -27 -27 27 -18 -9 18 9 18 9 -18 -9 -18 18 -9 9 18 -18 9 -9 -18 18 18 -18 -9 9 9 -9 -12 -6 -6 -3 12 6 6 3 -12 -6 12 6 -6 -3 6 3 -12 12 -6 6 -6 6 -3 3 -8 -4 -4 -2 -4 -2 -2 -1;
    18 -18 -18 18 -18 18 18 -18 9 9 -9 -9 -9 -9 9 9 12 -12 6 -6 -12 12 -6 6 12 -12 -12 12 6 -6 -6 6 6 6 3 3 -6 -6 -3 -3 6 6 -6 -6 3 3 -3 -3 8 -8 4 -4 4 -4 2 -2 4 4 2 2 2 2 1 1;
    -6 0 6 0 6 0 -6 0 0 0 0 0 0 0 0 0 -3 0 -3 0 3 0 3 0 -4 0 4 0 -2 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -2 0 -1 0 -1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 -6 0 6 0 6 0 -6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 -3 0 3 0 3 0 -4 0 4 0 -2 0 2 0 0 0 0 0 0 0 0 0 -2 0 -2 0 -1 0 -1 0;
    18 -18 -18 18 -18 18 18 -18 12 6 -12 -6 -12 -6 12 6 9 -9 9 -9 -9 9 -9 9 12 -12 -12 12 6 -6 -6 6 6 3 6 3 -6 -3 -6 -3 8 4 -8 -4 4 2 -4 -2 6 -6 6 -6 3 -3 3 -3 4 2 4 2 2 1 2 1;
    -12 12 12 -12 12 -12 -12 12 -6 -6 6 6 6 6 -6 -6 -6 6 -6 6 6 -6 6 -6 -8 8 8 -8 -4 4 4 -4 -3 -3 -3 -3 3 3 3 3 -4 -4 4 4 -2 -2 2 2 -4 4 -4 4 -2 2 -2 2 -2 -2 -2 -2 -1 -1 -1 -1;
    2 0 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 2 0 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -6 6 0 0 6 -6 0 0 -4 -2 0 0 4 2 0 0 0 0 0 0 0 0 0 0 -3 3 0 0 -3 3 0 0 0 0 0 0 0 0 0 0 -2 -1 0 0 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    4 -4 0 0 -4 4 0 0 2 2 0 0 -2 -2 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 2 -2 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6 6 0 0 6 -6 0 0 0 0 0 0 0 0 0 0 -4 -2 0 0 4 2 0 0 0 0 0 0 0 0 0 0 -3 3 0 0 -3 3 0 0 -2 -1 0 0 -2 -1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 -4 0 0 -4 4 0 0 0 0 0 0 0 0 0 0 2 2 0 0 -2 -2 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 2 -2 0 0 1 1 0 0 1 1 0 0;
    -6 0 6 0 6 0 -6 0 0 0 0 0 0 0 0 0 -4 0 -2 0 4 0 2 0 -3 0 3 0 -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 -2 0 -1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 -6 0 6 0 6 0 -6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4 0 -2 0 4 0 2 0 -3 0 3 0 -3 0 3 0 0 0 0 0 0 0 0 0 -2 0 -1 0 -2 0 -1 0;
    18 -18 -18 18 -18 18 18 -18 12 6 -12 -6 -12 -6 12 6 12 -12 6 -6 -12 12 -6 6 9 -9 -9 9 9 -9 -9 9 8 4 4 2 -8 -4 -4 -2 6 3 -6 -3 6 3 -6 -3 6 -6 3 -3 6 -6 3 -3 4 2 2 1 4 2 2 1;
    -12 12 12 -12 12 -12 -12 12 -6 -6 6 6 6 6 -6 -6 -8 8 -4 4 8 -8 4 -4 -6 6 6 -6 -6 6 6 -6 -4 -4 -2 -2 4 4 2 2 -3 -3 3 3 -3 -3 3 3 -4 4 -2 2 -4 4 -2 2 -2 -2 -1 -1 -2 -2 -1 -1;
    4 0 -4 0 -4 0 4 0 0 0 0 0 0 0 0 0 2 0 2 0 -2 0 -2 0 2 0 -2 0 2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 4 0 -4 0 -4 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 2 0 -2 0 -2 0 2 0 -2 0 2 0 -2 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0;
    -12 12 12 -12 12 -12 -12 12 -8 -4 8 4 8 4 -8 -4 -6 6 -6 6 6 -6 6 -6 -6 6 6 -6 -6 6 6 -6 -4 -2 -4 -2 4 2 4 2 -4 -2 4 2 -4 -2 4 2 -3 3 -3 3 -3 3 -3 3 -2 -1 -2 -1 -2 -1 -2 -1;
    8 -8 -8 8 -8 8 8 -8 4 4 -4 -4 -4 -4 4 4 4 -4 4 -4 -4 4 -4 4 4 -4 -4 4 4 -4 -4 4 2 2 2 2 -2 -2 -2 -2 2 2 -2 -2 2 2 -2 -2 2 -2 2 -2 2 -2 2 -2 1 1 1 1 1 1 1 1]
)

const GAS_NAMES = Dict{String,String}(
    "Ar"    => "argon",
    "CH4"   => "methane",
    "CH4O"  => "methanol",
    "CH3OH" => "methanol",
    "CO"    => "carbon monoxide",
    "CO2"   => "carbon dioxide",
    "H2"    => "hydrogen",
    "H4Si"  => "silane",
    "SiH4"  => "silane",
    "H2O"   => "water",
    "He"    => "helium", # missing for PCSAFT
    "Kr"    => "krypton",
    "N2"    => "nitrogen",
    "NO"    => "nitrogen monoxide", # missing for PCSAFT
    "NO2"   => "nitrogen dioxide", # missing for PCSAFT
    "Ne"    => "neon", # missing for PCSAFT
    "O2"    => "oxygen",
    "O2S"   => "sulfur dioxide",
    "SO2"   => "sulfur dioxide",
    "Xe"    => "xenon",
)

const GERG2008_nameset = Set(["methane", "nitrogen", "carbon dioxide", "ethane", "propane", "butane", "isobutane", "pentane", "isopentane", "hexane", "heptane", "octane", "nonane", "decane", "hydrogen", "oxygen", "carbon monoxide", "water", "hydrogen sulfide", "helium", "argon"])
const 𝒩ₐ = 6.02214076e23u"mol^-1"

# color scheme
const atom_info = [
    (:H, "#FFFFFF", "#FFFFFF"),
    (:He, "#D9FFFF", "#FFC0CB"),
    (:Li, "#CC80FF", "#B22222"),
    (:Be, "#C2FF00", "#FF1493"),
    (:B, "#FFB5B5", "#00FF00"),
    (:C, "#909090", "#C8C8C8"),
    (:N, "#3050F8", "#8F8FFF"),
    (:O, "#FF0D0D", "#F00000"),
    (:F, "#90E050", "#DAA520"),
    (:Ne, "#B3E3F5", "#FF1493"),
    (:Na, "#AB5CF2", "#0000FF"),
    (:Mg, "#8AFF00", "#228B22"),
    (:Al, "#BFA6A6", "#808090"),
    (:Si, "#F0C8A0", "#DAA520"),
    (:P, "#FF8000", "#FFA500"),
    (:S, "#FFFF30", "#FFC832"),
    (:Cl, "#1FF01F", "#00FF00"),
    (:Ar, "#80D1E3", "#FF1493"),
    (:K, "#8F40D4", "#FF1493"),
    (:Ca, "#3DFF00", "#808090"),
    (:Sc, "#E6E6E6", "#FF1493"),
    (:Ti, "#BFC2C7", "#808090"),
    (:V, "#A6A6AB", "#FF1493"),
    (:Cr, "#8A99C7", "#808090"),
    (:Mn, "#9C7AC7", "#808090"),
    (:Fe, "#E06633", "#FFA500"),
    (:Co, "#F090A0", "#FF1493"),
    (:Ni, "#50D050", "#A52A2A"),
    (:Cu, "#C88033", "#A52A2A"),
    (:Zn, "#7D80B0", "#A52A2A"),
    (:Ga, "#C28F8F", "#FF1493"),
    (:Ge, "#668F8F", "#FF1493"),
    (:As, "#BD80E3", "#FF1493"),
    (:Se, "#FFA100", "#FF1493"),
    (:Br, "#A62929", "#A52A2A"),
    (:Kr, "#5CB8D1", "#FF1493"),
    (:Rb, "#702EB0", "#FF1493"),
    (:Sr, "#00FF00", "#FF1493"),
    (:Y, "#94FFFF", "#FF1493"),
    (:Zr, "#94E0E0", "#FF1493"),
    (:Nb, "#73C2C9", "#FF1493"),
    (:Mo, "#54B5B5", "#FF1493"),
    (:Tc, "#3B9E9E", "#FF1493"),
    (:Ru, "#248F8F", "#FF1493"),
    (:Rh, "#0A7D8C", "#FF1493"),
    (:Pd, "#006985", "#FF1493"),
    (:Ag, "#C0C0C0", "#808090"),
    (:Cd, "#FFD98F", "#FF1493"),
    (:In, "#A67573", "#FF1493"),
    (:Sn, "#668080", "#FF1493"),
    (:Sb, "#9E63B5", "#FF1493"),
    (:Te, "#D47A00", "#FF1493"),
    (:I, "#940094", "#A020F0"),
    (:Xe, "#429EB0", "#FF1493"),
    (:Cs, "#57178F", "#FF1493"),
    (:Ba, "#00C900", "#FFA500"),
    (:La, "#70D4FF", "#FF1493"),
    (:Ce, "#FFFFC7", "#FF1493"),
    (:Pr, "#D9FFC7", "#FF1493"),
    (:Nd, "#C7FFC7", "#FF1493"),
    (:Pm, "#A3FFC7", "#FF1493"),
    (:Sm, "#8FFFC7", "#FF1493"),
    (:Eu, "#61FFC7", "#FF1493"),
    (:Gd, "#45FFC7", "#FF1493"),
    (:Tb, "#30FFC7", "#FF1493"),
    (:Dy, "#1FFFC7", "#FF1493"),
    (:Ho, "#00FF9C", "#FF1493"),
    (:Er, "#00E675", "#FF1493"),
    (:Tm, "#00D452", "#FF1493"),
    (:Yb, "#00BF38", "#FF1493"),
    (:Lu, "#00AB24", "#FF1493"),
    (:Hf, "#4DC2FF", "#FF1493"),
    (:Ta, "#4DA6FF", "#FF1493"),
    (:W, "#2194D6", "#FF1493"),
    (:Re, "#267DAB", "#FF1493"),
    (:Os, "#266696", "#FF1493"),
    (:Ir, "#175487", "#FF1493"),
    (:Pt, "#D0D0E0", "#FF1493"),
    (:Au, "#FFD123", "#DAA520"),
    (:Hg, "#B8B8D0", "#FF1493"),
    (:Tl, "#A6544D", "#FF1493"),
    (:Pb, "#575961", "#FF1493"),
    (:Bi, "#9E4FB5", "#FF1493"),
    (:Po, "#AB5C00", "#FF1493"),
    (:At, "#754F45", "#FF1493"),
    (:Rn, "#428296", "#FF1493"),
    (:Fr, "#420066", "#FF1493"),
    (:Ra, "#007D00", "#FF1493"),
    (:Ac, "#70ABFA", "#FF1493"),
    (:Th, "#00BAFF", "#FF1493"),
    (:Pa, "#00A1FF", "#FF1493"),
    (:U, "#008FFF", "#FF1493"),
    (:Np, "#0080FF", "#FF1493"),
    (:Pu, "#006BFF", "#FF1493"),
    (:Am, "#545CF2", "#FF1493"),
    (:Cm, "#785CE3", "#FF1493"),
    (:Bk, "#8A4FE3", "#FF1493"),
    (:Cf, "#A136D4", "#FF1493"),
    (:Es, "#B31FD4", "#FF1493"),
    (:Fm, "#B31FBA", "#FF1493"),
    (:Md, "#B30DA6", "#FF1493"),
    (:No, "#BD0D87", "#FF1493"),
    (:Lr, "#C70066", "#FF1493"),
    (:Rf, "#CC0059", "#FF1493"),
    (:Db, "#D1004F", "#FF1493"),
    (:Sg, "#D90045", "#FF1493"),
    (:Bh, "#E00038", "#FF1493"),
    (:Hs, "#E6002E", "#FF1493"),
    (:Mt, "#EB0026", "#FF1493"),
]


# other utils
nint(x) = floor(Int, ifelse(x>=0.0, x+0.5, x-0.5))
