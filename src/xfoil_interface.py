import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt


class Xfoil:

    def __init__(self, xfoil_location, N_iter, Re, M, alfas):
        self.xfoil_dir = xfoil_location
        self.N = N_iter
        self.Re = Re
        self.M = M
        self.alfas = alfas

    def issueCmd(self, ps, cmd):
        """
        issue a command into the shell
        :param cmd: command to execture (string)
        """

        ps.stdin.write((cmd + '\n').encode('utf-8'))

    def get_polar(self, airfoil_dat, polar_txt):
        """
        Generate an airfoil polar by running an xfoil analysis.
        :param airfoil_dat: .dat file containing airfoil coordinates
        :param polar_txt: name of .txt file to store xfoil output
        """
        ps = sp.Popen([self.xfoil_dir], stdin=sp.PIPE, stderr=sp.PIPE, stdout=sp.PIPE)

        # if os.path.exists(polar_txt):
        #     os.remove(polar_txt)

        self.issueCmd(ps, 'load ' + airfoil_dat)
        self.issueCmd(ps, 'oper')
        self.issueCmd(ps, 'visc ' + str(self.Re))
        self.issueCmd(ps, 'M ' + str(self.M))
        self.issueCmd(ps, 'iter ' + str(self.N))
        self.issueCmd(ps, 'pacc')
        self.issueCmd(ps, polar_txt)
        self.issueCmd(ps, ' ')
        a1, a2, inc = [str(x) for x in self.alfas]
        self.issueCmd(ps, 'aseq ' + ' ' + a1 + ' ' + a2 + ' ' + inc)

        ps.communicate('quit'.encode('utf-8'))

    def get_opt_params(self, polar_txt):
        """
        :param polar_txt: file containing the output of xfoil polar analysis.
        :return: maximum CL3 / CD2 ratio for the set of angles of attack.
        """
        with open(polar_txt) as file:
            polar_data = np.array([np.array([float(x) for x in line.split()]) for line in file.readlines()[12:]])
            alfa = polar_data[:,0]
            Cl = polar_data[:,1]
            Cd = polar_data[:,2]

            idx = np.argmax(Cl*Cl*Cl / Cd / Cd)
            ClCd = Cl[idx] / Cd[idx]
            Cl3Cd2 = (Cl[idx])**3 / (Cd[idx])**2

            stall_idx = np.argmax(Cl)
            alfa_range = alfa[stall_idx] - alfa[idx]

        return ClCd, Cl3Cd2, alfa_range

    def plot_polar(self, polar_txt):
        """
        Plot Cl/alpha, Cd/alpha and Cl/Cd and Cl^3/Cd^2 from a certain analysis.
        :param polar_txt: .txt file where xfoil polar data is stored.
        """
        with open(polar_txt) as file:
            data = np.array([np.array([float(x) for x in line.split()]) for line in file.readlines()[12:]])
            alfa = data[:,0]
            Cl = data[:,1]
            Cd = data[:,2]
            Cl3Cd2 = Cl**3 / Cd**2

        fig, axs = plt.subplots(2,2)

        axs[0, 0].plot(alfa, Cl, marker='x', color='limegreen')
        axs[0, 0].set(xlabel=r'$\alpha$ [-]', ylabel='$C_{l}$ [-]')

        axs[0, 1].plot(alfa, Cd, marker='x', color='mediumblue')
        axs[0, 1].set(xlabel=r'$\alpha$ [-]', ylabel='$C_{d}$ [-]')

        axs[1, 0].plot(Cd, Cl, marker='x', color='orangered')
        axs[1, 0].set(xlabel=r'$C_{d}$ [-]', ylabel='$C_{l}$ [-]')

        axs[1, 1].plot(alfa, Cl3Cd2, marker='x', color='black')
        axs[1, 1].set(xlabel=r'$\alpha$ [-]', ylabel='$C_{l}^{3}/C_{d}^{2}$ [-]')

        plt.show()


if __name__ == "__main__":
    plt.style.use('classic')
    plt.rcParams.update({'font.size': 20})

    airfoil_file = 'airfoils/nlf1015.dat'
    polar_file = 'airfoils/polar.dat'
    xfoil_dir = '/bin/xfoil'

    Re = 8.0955e5
    M = 0.5
    N = 400
    alfas = [0, 13, 0.5]

    xfoil = Xfoil(xfoil_dir, N, Re, M, alfas)
    xfoil.get_polar(airfoil_file, polar_file)
    xfoil.plot_polar(polar_file)
    #params = xfoil.get_opt_params(polar_file)

