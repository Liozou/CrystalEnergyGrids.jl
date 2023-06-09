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

# Van der Waals state equation coefficients (from the Handbook (?))
const VDW_COEFF = Dict{String,Tuple{Float64,Float64}}(
    "H2"    => (0.2476, 0.02661),
    "CH4"   => (2.283, 0.04278),
    "CH3OH" => (9.649, 0.06702),
    "CO"    => (1.505, 0.03985),
    "CO2"   => (3.640, 0.04267),
    "H2O"   => (5.536, 0.03049),
    "N2"    => (1.408, 0.03913),
    "NO"    => (1.358, 0.02789),
    "NO2"   => (5.354, 0.04424),
    "O2"    => (1.378, 0.03183),
    "SiH4"  => (4.377, 0.05786),
    "SO2"   => (6.803, 0.05636),
    "He"    => (0.03457, 0.0237),
    "Ne"    => (0.2135, 0.01709),
    "Ar"    => (1.363, 0.03219),
    "Kr"    => (2.349, 0.03978),
    "Xe"    => (4.25, 0.05105),
)

# other utils
nint(x) = floor(Int, ifelse(x>=0.0, x+0.5, x-0.5))
