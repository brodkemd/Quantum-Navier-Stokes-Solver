
import numpy as np

def Calc_dffdt(dFdt, dJdt, Del_x, d, Tot_Int_Pts):
    #CALC_DFFDT evaluates first time derivative of ODE driver function
    #   Calc_dffdt evaluates first time derivative of ODE driver function at
    #       all interior grid-points.
    #
    #   INPUTS:
    #           dFdt = d x Tot_X_Pts array storing first time derivative of
    #                       flow flux at all grid-points
    #           dJdt = 1 x Tot_Int_Pts array storing first time derivative of
    #                       flow source term at all interior grid-points
    #           Del_x = distance between grid-points
    #           d = number of primary flow variables at a grid-point
    #           Tot_Int_pts = number of interior grid-points
    #
    #   OUTPUT:
    #           dffdt_vals = d x Tot_Int_Pts array storing first time
    #                           derivative of ODE driver function

    dffdt_vals = np.zeros((d,Tot_Int_Pts))

    TotIPtp1 = Tot_Int_Pts + 1

    fac1 = (-1)/(2*Del_x)

    # evaluate dffdt_vals at interior points (2 <= ll <= TotIPtp1)

    for ll in range(1,TotIPtp1):
        IPLabel = ll - 1

        dffdt_vals[0, IPLabel] = fac1*(dFdt[0, ll+1] - dFdt[0, ll-1])

        dffdt_vals[1, IPLabel] = fac1*(dFdt[1, ll+1] - dFdt[1, ll-1]-dJdt[IPLabel])

        dffdt_vals[2, IPLabel] = fac1*(dFdt[2, ll+1] - dFdt[2, ll-1])

    return dffdt_vals
