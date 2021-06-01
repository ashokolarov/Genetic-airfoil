from Airfoil import Airfoil
from xfoil_interface import Xfoil
import numpy as np
import random
import matplotlib.pyplot as plt
import multiprocessing as mp
from tabulate import tabulate
import os
import shutil


def cost_function(parsec):
    """
    Obtain cost of airfoil from its parsec parameters.
    Cost function = 0.5*Cl3Cd2 + 0.3*ClCd + 0.2*alfa_range
    :param parsec: Parsec parameters of airfoil
    :return: cost
    """
    name = int(mp.current_process()._identity[0])
    optim_file = f'polars/optimized{name}.dat'
    polar_file = f'polars/polar{name}.dat'

    _foil = Airfoil.from_parsec(parsec)
    _foil.save_as_dat(optim_file, 150)
    xfoil.get_polar(optim_file, polar_file)

    _clcd, _cl3cd2, _alfa_range = xfoil.get_opt_params(polar_file)

    _cost = 0.5*(_cl3cd2/cl3cd2_0) + 0.3*(_clcd/clcd_0) + 0.2*(_alfa_range/alfa_range_0)
    return _cost


def crossover(_parents, prob_cross):
    """
    Generate cross over children from two parents depended on some crossover probability.
    :param _parents: list containing the two parents.
    :param prob_cross: probability of crossover
    :return: two produced children
    """
    c1, c2 = _parents[0].copy(), _parents[1].copy()
    if random.random() < prob_cross:
        pt = random.randint(1, param_size - 2)
        c1[:pt] = _parents[0][:pt]
        c1[pt:] = _parents[1][pt:]
        c2[:pt] = _parents[1][:pt]
        c2[pt:] = _parents[0][pt:]
    return c1, c2


def mutate(x, prob_mutate):
    """
    Perform mutation on a gene dependent on some probability.
    :param x: parsec parameters of the child
    :param prob_mutate: probability of mutation
    :return: mutated parsec parameters.
    """
    for k in range(param_size):
        if random.random() < prob_mutate:
            x[k] = x0[k] * (1 + random.choice([-MAX_CHANGE, MAX_CHANGE]))


if __name__ == "__main__":
    xfoil_dir = '/bin/xfoil'
    original_airfoil = 'airfoils/nlf1015.dat'
    original_polar = 'airfoils/polar.dat'

    if not os.path.exists('polars'):
        os.mkdir('polars')

    # Generate initial airfoil to obtain param and cost
    foil = Airfoil.from_dat(original_airfoil)
    x0 = foil.parsec_params

    # Set xfoil parameters
    Re = 8.0955e5
    M = 0.5
    N = 1000
    alfas = [0, 12, 0.5]

    # Generate xfoil instance
    xfoil = Xfoil(xfoil_dir, N, Re, M, alfas)

    # Obtain cost parameters for initial airfoil
    clcd_0, cl3cd2_0, alfa_range_0 = xfoil.get_opt_params(original_polar)  # Aero parameters of the initial airfoil

    MAX_CHANGE = 0.1       # Maximum change in parameters

    # GENETIC
    param_size = len(x0)   # Number of parameters to optimize
    init_pop_size = 60     # Size of initial population
    mating_pool_size = 20  # Number of parents to mate
    p_cross = 0.8          # Probability of crossover
    p_mutate = 0.025       # Probability of mutation
    elitism = 0.05         # Proportion of parents to carry into next generation
    iter_num = 15          # Number of iterations

    # Generate initial population
    population = np.zeros((init_pop_size, param_size))

    for i in range(init_pop_size):
        sol = np.zeros(param_size)
        for j in range(param_size):
            sol[j] = x0[j] * (1 + np.random.uniform(-MAX_CHANGE, MAX_CHANGE))
        population[i] = sol

    best_x, best_val = x0, 1
    cost_progression = [[1]]

    # Start genetic algorithm
    num_processes = mp.cpu_count()
    pool = mp.Pool(processes=num_processes)
    for cur_iter in range(iter_num):
        # Evaluate pool
        evals = list(pool.map(cost_function, population))
        selected = [x for x, _ in sorted(zip(population, evals), key=lambda pair: -pair[1])][:mating_pool_size]
        cost_progression.append(evals)

        # Check if a new best solution is found
        if max(evals) > best_val:
            best_x = selected[0]
            best_val = max(evals)

        # Create list of childrean
        children = list()

        # Apply elitism
        selection_size = len(selected)
        parents_left = int(elitism * selection_size)
        for i in range(parents_left):
            children.append(selected[i])

        # Choose parents with weighted random func based on their fitness
        w = [selection_size/(i+selection_size) for i in range(selection_size)]

        # Perform crossover and mutation
        for i in range(0, len(population), 2):
            parents = random.choices(selected, weights=w, k=2)

            for c in crossover(parents, p_cross):
                mutate(c, p_mutate)
                children.append(c)

        population = children
        print(f'Iteration number: {cur_iter}')
        print(f'Current best cost: {best_val}\n')

    optimized = Airfoil.from_parsec(best_x)
    optimized.save_as_dat("airfoils/optimized.dat", 150)
    xfoil.get_polar("airfoils/optimized.dat", "airfoils/optimized_polar.dat")

    clcd, cl3cd2, alfa_range = xfoil.get_opt_params('airfoils/optimized_polar.dat')

    print('Genetic algorithm finished')
    table = tabulate([['Parameter', 'Initial airfoil', 'Optimized airfoil'],
             ['ClCd', clcd_0, clcd],
             ['Cl3Cd2', cl3cd2_0, cl3cd2],
             ['Angle of attack range', alfa_range_0, alfa_range]], headers='firstrow')
    print(table)

    initial = Airfoil.from_dat('airfoils/nlf1015.dat')
    x = Airfoil.grid_chebychev(100)

    plt.style.use('ggplot')
    fig, ax = plt.subplots(1)

    ax.plot(x, [optimized.zup(i) for i in x], color='orange', label='optimized')
    ax.plot(x, [optimized.zlo(i) for i in x], color='orange')
    ax.plot(x, [initial.zup(i) for i in x], color='blue', label='initial')
    ax.plot(x, [initial.zlo(i) for i in x], color='blue')

    ax.set_aspect('equal')
    fig.tight_layout()
    plt.xlabel('Normalized chord [-]')
    plt.ylabel('Normalized thickness [-]')
    plt.grid(True)
    plt.legend()
    plt.show()

    with open('cost_progression.txt', 'w') as f:
        for i, elem in enumerate(cost_progression):
            f.write(f"{i} : ")
            for x in elem:
                f.write(f"{x}, ")
            f.write('\n')

    for i, cost in enumerate(cost_progression):
        x = [i for x in cost]
        y = cost
        plt.scatter(x, y, marker='o', c=y, cmap='RdYlGn')
    max_vals = [max(x) for x in cost_progression]
    plt.plot(max_vals, linestyle='--', color='blue')
    plt.xlabel('Iteration number [-]')
    plt.ylabel('Cost [-]')
    plt.show()

    shutil.rmtree('polars')
