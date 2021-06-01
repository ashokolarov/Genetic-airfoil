import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize


class Airfoil:

    def __init__(self, parsec_params, zup, zlo):
        """
        Airfoil class to translate between coordinates and PARSEC elements.

        :param parsec_params: array holding parsec parameters
        :param zup: function z(x) returning z-coordinate of upper surface.
        :param zlo: function z(x) returning z-coordinate of lower surface.
        """
        self.parsec_params = parsec_params
        self.zup = zup
        self.zlo = zlo
        self.labels = ['rle', 'Xup', 'Zup', 'Zxxup', 'Xlo', 'Zlo',
                       'Zxxlo', 'Zte', 'DZte', 'alte', 'bete']

    @classmethod
    def from_parsec(cls, *args):
        """
        Create a class instance from number of nodes and PARSEC elements.
        :param args: array holding parsec parameters
        """
        Aup, Alo = Airfoil._get_poly_coeff(*args)

        def zup(x):
            _zup = 0
            for i in range(6):
                _zup += Aup[i] * x **(i+0.5)
            return _zup

        def zlo(x):
            _zlo = 0
            for i in range(6):
                _zlo += Alo[i] * x **(i+0.5)
            return _zlo

        return Airfoil(*args, zup, zlo)

    @staticmethod
    def _get_poly_coeff(args):
        """
        Find polynomial coefficients of interpolant from PARSEC elements.

        :param args: array of Parsec elements.
        """
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11 = args

        p10, p11 = np.deg2rad([p10, p11])

        coeff = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]

        # Upper surface 
        Cup = np.ones((6,6), dtype=np.float64)

        for i in range(6):
            Cup[1,i] = p2 ** coeff[i]
            Cup[2,i] = coeff[i]
            Cup[3,i] = coeff[i] * p2 **(coeff[i] - 1)
        
        Cup[4,0] = -0.25 * p2 **(-1.5)
        Cup[4,1] =  0.75 * p2 **(-0.5)
        Cup[4,2] =  3.75 * p2 **( 0.5)
        Cup[4,3] =  8.75 * p2 **( 1.5)
        Cup[4,4] = 15.75 * p2 **( 2.5)
        Cup[4,5] = 24.75 * p2 **( 3.5)
        
        Cup[5,1:] = 0

        Bup = np.zeros(6)

        Bup[0] = p8 + p9/2
        Bup[1] = p3
        Bup[2] = np.tan(p10 - p11/2)
        Bup[4] = p4
        Bup[5] = np.sqrt(2*p1)

        Aup = np.linalg.solve(Cup, Bup)
        
        # Lower surface
        Clo = np.ones((6,6), dtype=np.float64)

        for i in range(6):
            Clo[1,i] = p5 ** coeff[i]
            Clo[2,i] = coeff[i]
            Clo[3,i] = coeff[i] * p5 **(coeff[i] - 1)

        Clo[4,0] = -0.25 * p5 **(-1.5)
        Clo[4,1] =  0.75 * p5 **(-0.5)
        Clo[4,2] =  3.75 * p5 **( 0.5)
        Clo[4,3] =  8.75 * p5 **( 1.5)
        Clo[4,4] = 15.75 * p5 **( 2.5)
        Clo[4,5] = 24.75 * p5 **( 3.5)

        Clo[5,1:] = 0

        Blo = np.zeros(6)
        
        Blo[0] = p8 - p9/2
        Blo[1] = p6
        Blo[2] = np.tan(p10 + p11/2)
        Blo[4] = p7
        Blo[5] = -np.sqrt(2*p1)

        Alo = np.linalg.solve(Clo, Blo)

        return Aup, Alo

    def print_parsec(self):
        """
        Print PARSEC elements of the airfoil
        """
        for i in range(len(self.labels)):
            print(f"{self.labels[i]} : {self.parsec_params[i]:.3f}")

    @classmethod
    def from_dat(cls, dat_file, plot=False):
        """
        Obtain parsec parameters from an airfoil .dat file by solving an optimization problem.

        :param dat_file: .dat file holding coordinates of airfoil.
        """
        with open(dat_file, 'r') as points:
            x, y = np.loadtxt(points, unpack=True, skiprows=1)
        xy = np.array(list(zip(x,y)))

        def cost(x0):
            foil = Airfoil.from_parsec(x0)

            y_fit = np.zeros(len(x))
            
            for i in range(len(x)):
                if i<32:
                    y_fit[i] = foil.zup(x[i])
                else:
                    y_fit[i] = foil.zlo(x[i])

            err = np.sum((y - y_fit) ** 2)
            return err

        # Set up optimization problem
        x0 = [0.0083, 0.423, 0.0587, -0.347, 0.358, -0.032, 0.417, 0, 0, 10.03, 5.64]
        bnd = ((0, None), (0, None), (0, None), (None, None), (0, None), (None, None), (None, None), (-0.00001, 0.00001),
               (-0.00001, 0.00001), (-45, 45), (-45, 45))
        res = minimize(cost, x0, method='SLSQP', bounds=bnd, tol=1e-16).x
        res[7] = 0
        res[8] = 0

        if plot:
            foil = Airfoil.from_parsec(res)

            y_fit = np.zeros(len(x))

            for i in range(len(x)):
                if i < 32:
                    y_fit[i] = foil.zup(x[i])
                else:
                    y_fit[i] = foil.zlo(x[i])

            fig, ax = plt.subplots(1)
            ax.plot(x, y, label='True')
            ax.plot(x, y_fit, label='Fitted')
            ax.set_aspect('equal')
            fig.tight_layout()
            plt.grid(True)
            plt.xlabel('Normalized chord [-]')
            plt.ylabel('Normalized thickness [-]')
            plt.legend()
            plt.show()
            
        return cls.from_parsec(res)

    @staticmethod
    def grid_chebychev(n):
        """
        Generate Chebychev grid over [0,1] interval.
        :param n: number of grid points.
        """
        a = 0
        b = 1
        xi = np.cos( np.pi/(2*(n+1)) * (2*np.linspace(1, n+1, n+1) - 1) )
        return ((-xi + 1)/2) * (b-a) + a

    def save_as_dat(self, name, N, plot=False):
        """
        Save generated airfoil coordinates to a .dat file using chebychev nodes.
        :param name: name of file to save to
        :param N: number of total nodes to use
        """
        x = self.grid_chebychev(N)
        zup = [self.zup(i) for i in x]
        zlo = [self.zlo(i) for i in x]

        if plot:
            plt.plot(x, zup, marker='x')
            plt.plot(x, zlo, marker='x')
            plt.show()

        data = np.array(list(zip(x, zup)))
        data = np.flip(data,0)

        data_lo = np.array(list(zip(x, zlo)))
        data = np.vstack((data, data_lo))

        if '.dat' not in name:
            name += '.dat'

        header = "Airfoil"
        np.savetxt(name, data, header=header, comments="")

    def plot(self):
        """
        Plot airfoil geometry.
        """
        N = 50
        x = self.grid_chebychev(N)
        zup = [self.zup(i) for i in x]
        zlo = [self.zlo(i) for i in x]

        fig, ax = plt.subplots(1)
        
        ax.plot(x, zup, marker='x')
        ax.plot(x, zlo, marker='o')
        
        ax.set_aspect('equal')
        fig.tight_layout()
        plt.xlabel('Normalized chord [-]')
        plt.ylabel('Normalized thickness [-]')
        plt.grid(True)
        
        plt.show()
        

if __name__ == "__main__":
    dat_file = "nlf1015.dat"
    foil = Airfoil.from_dat(dat_file)
    foil.save_as_dat('original.dat', 150, True)
    
    #save_as = 'bla.dat'
    #N = 150
    #foil.save_as_dat(save_as, N)

