

# definitions from RASPA

const ANGSTROM = 1e-10
const NANO_SECOND = 1e-9
const PICO_SECOND = 1e-12
const ATOMIC_MASS_UNIT = 1.6605402e-27
const MOLAR_GAS_CONSTANT = 8.314464919
const AVOGADRO_CONSTANT = 6.0221419947e23
const ENERGY_CONVERSION_FACTOR = ATOMIC_MASS_UNIT * ANGSTROM^2/PICO_SECOND^2
const ENERGY_TO_KELVIN = ENERGY_CONVERSION_FACTOR * AVOGADRO_CONSTANT / MOLAR_GAS_CONSTANT
const ELECTRONIC_CHARGE_UNIT = 1.60217733e-19
const ELECTRIC_CONSTANT = 8.8541878176e-12
const DielectricConstantOfTheMedium = 1.0
const COULOMBIC_CONVERSION_FACTOR = ELECTRONIC_CHARGE_UNIT^2/(4π*ELECTRIC_CONSTANT*ANGSTROM*ENERGY_CONVERSION_FACTOR*DielectricConstantOfTheMedium)


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

# definitions from Unitful
const ANG_UNIT = 1.0AtomsBase.Unitful.Å
const CHARGE_UNIT = 1.0AtomsBase.UnitfulAtomic.e_au
const ATOMMASS_UNIT = 1.0AtomsBase.Unitful.u


# other utils
nint(x) = floor(Int, ifelse(x>=0.0, x+0.5, x-0.5))
