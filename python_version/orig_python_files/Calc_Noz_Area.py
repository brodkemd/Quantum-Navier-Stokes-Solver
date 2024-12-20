import numpy as np

def Calc_Noz_Area( x ):
    #CALC_NOZ_AREA determines nozzle area at grid-points x
    #   Calc_Noz_Area evaluates the nozzle area at grid-points x.
    #       
    #   INPUT:
    #           x = array with spatial grid-point locations 
    #
    #   OUTPUT:
    #           A = 1 x Tot_Pts array with nozzle area at each grid-point x
    #

    # determine number of grid-points

    Tot_Pts = len(x)

    A = np.zeros((Tot_Pts))

    # Nozzle throat located at x = 1.5 in formula below
    for i in range(Tot_Pts):
        A[i] = 1 + 2.2*( x[i] - 1.5 )**(2)

    return A