import numpy as np
from acados_settings import *




def main():
    """ Main program """
    Tsim = 20
    Tf = 1
    N = 50
    Nsim = N/Tf*Tsim

    constraints, model, acados_solver = acados_settings(Tf, N)

    return 0

if __name__ == "__main__":
    main()
