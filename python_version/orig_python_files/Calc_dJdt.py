
import numpy as np

def Calc_dJdt(U, ff_vals, A, Gamma, d, Tot_Int_Pts):
    #CALC_DJDT evaluates first time derivative of source term at interior pts
    #   Calc_dJdt evaluates first time derivative of flow source term for
    #       1D compressible, inviscid flow through nozzle at all interior
    #       grid-points.
    #
    #   INPUTS:
    #           U = d x Tot_X_Pts array storing primary flow variables at all
    #                               grid-points
    #           ff_vals = d x Tot_Int_Pts array storing ODE driver function
    #                               values at all interior grid-points
    #           A = 1 x Tot_X_Pts array with nozzle area at all grid-points
    #           Gamma = ratio of specific heats
    #           d = number of primary flow variables at each grid-point
    #           Tot_Int_Pts = number of interior grid-points
    #
    #   OUTPUT:
    #           dJdt = 1 x Tot_Int_Pts array storing first time derivative of
    #                                   flow source term at all interior
    #                                   grid-points.

    dJdt = np.zeros(Tot_Int_Pts)

    TotIPtp1 = Tot_Int_Pts + 1

    fac1 = (Gamma - 1)/Gamma
    fac2 = Gamma/2

    # evaluate dJdt by looping over interior grid-points 
    #           ( 2 <= ll < = TotIPtp1 )

    for ll in range(TotIPtp1):
        IPLabel = ll - 1

        fac3 = U[1, ll]/U[0, ll]
        fac4 = np.log(A[ll+1]/A[ll-1])

        dJdt[IPLabel] = fac4*fac1*(ff_vals[2, IPLabel]-fac2*fac3*(2*ff_vals[1, IPLabel]-fac3*ff_vals[0, IPLabel] ) )

    return dJdt