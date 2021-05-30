from Airfoil import Airfoil
from xfoil_interface import Xfoil


def cost_function(parsec):
    foil = Airfoil.from_parsec(parsec)
    foil.save_as_dat(optim_file, 150)
    xfoil.get_polar(optim_file, polar_file)
    cl3cd2 = xfoil.get_CL3CD2(polar_file)
    return cl3cd2

if __name__ == "__main__":
    xfoil_dir = '/bin/xfoil'
    original_airfoil = 'airfoils/nlf1015.dat'

    foil = Airfoil.from_dat(original_airfoil)
    x0 = foil.parsec_params

    optim_file = 'airfoils/optimized.dat'
    polar_file = 'airfoils/polar.dat'

    Re = 8.0955e5
    M = 0.5
    N = 400
    alfas = [0, 10, 0.5]

    xfoil = Xfoil(xfoil_dir, N, Re, M, alfas)
    cl3cd2 = cost_function(x0)

    # Rle Xup Zup Zxxup Xlo Zlo Zxxlo Zte DZte alte bte

    MAX_CHANGE = 0.05
    bnds = [(x * (1 - MAX_CHANGE), x * (1 + MAX_CHANGE)) for x in x0]
