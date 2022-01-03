import numpy as np

def Calc_dFdt(U,ff_vals,ff_Bvals,Gamma,d,Tot_X_Pts):
    #CALC_DFDT evaluates first time derivative of flow fluxes at grid-points
    #   Calc_dFdt evaluates first time derivative of flow fluxes at
    #   all grid-points for 1D compressible, inviscid flow through a nozzle.
    #
    #   INPUTS:
    #           U = d x Tot_X_Pts array with primary flow variables at all
    #                   grid-points
    #           ff_vals = d x Tot_Int_Pts array with value of ODE driver
    #                       function at all interior grid-points
    #           ff_Bvals = d x 2 array with values of dU/dt at boundary points
    #           Gamma = ratio of specific heats
    #           d = number of primary flow variables at a grid-point
    #           Tot_X_Pts = number of grid-points
    #
    #   OUTPUT:
    #           dFdt = d x Tot_X_Pts array storing first time derivative of
    #                       flow fluxes at all grid-points

    dFdt = np.zeros((d, Tot_X_Pts))

    TotXPtm1 = Tot_X_Pts - 1

    fac4 = Gamma -1
    fac5 = Gamma/2
    fac6 = fac4*fac5
    fac7 = (3-Gamma)/2
    fac8 = fac4/Gamma

    # evaluate dFdt at boundaries:

    for ll in range(2):
        if ll == 0: # at nozzle entrance
            fac1 = U[1, 0]/U[0, 0]
            fac2 = U[2, 0]/U[0, 0]
            fac3 = fac1*fac2
            fac9 = (fac1)**(2)

            dFdt[0, ll] = ff_Bvals[2, ll]

            dFdt[1, ll] = fac7*fac1*(2*ff_Bvals[1, ll]-fac1*ff_Bvals[0, ll]) + fac8*ff_Bvals[2, ll]

            dFdt[2, ll] = Gamma*(fac1*ff_Bvals[2, ll]+fac2*ff_Bvals[1, ll]-fac3*ff_Bvals[0, ll])-fac6*fac9*(3*ff_Bvals[1, ll])-2*fac1*ff_Bvals[0, ll]
        elif ll == 1: # at nozzle exit
            fac1 = U[1, Tot_X_Pts - 1]/U[0, Tot_X_Pts - 1]
            fac2 = U[2, Tot_X_Pts - 1]/U[0, Tot_X_Pts - 1]
            fac3 = fac1*fac2
            fac9 = (fac1)**(2)

            dFdt[0, Tot_X_Pts- 1] = ff_Bvals[1, ll]

            dFdt[1, Tot_X_Pts- 1]= fac7*fac1*(2*ff_Bvals[1, ll]-fac1*ff_Bvals[0, ll]) + fac8*ff_Bvals[2, ll]

            dFdt[2, Tot_X_Pts- 1] = Gamma*(fac1*ff_Bvals[2, ll] + fac2*ff_Bvals[2, ll]- fac3*ff_Bvals[0, ll])-fac6*fac9*(3*ff_Bvals[1, ll]-2*fac1*ff_Bvals[0, ll])
        else:
            print('Unknown switch case label: ', ll)

    # evaluate dFdt at interior grid-points

    for ll in range(1, TotXPtm1 + 1):
        fac1 = U[1, ll]/U[0, ll]
        fac2 = U[2, ll]/U[0, ll]
        fac3 = fac1*fac2
        fac9 = (fac1)**(2)

        IPLabel = ll -1

        dFdt[0, ll] = ff_vals[1, IPLabel- 1]

        dFdt[1, ll] = fac7*fac1*(2*ff_vals[1, IPLabel- 1]-fac1*ff_vals[0, IPLabel- 1])+ fac8*ff_vals[2, IPLabel- 1]

        dFdt[2, ll] = Gamma*(fac1*ff_vals[2, IPLabel- 1] + fac2*ff_vals[1, IPLabel- 1] - fac3*ff_vals[1, IPLabel- 1])-fac6*fac9*(3*ff_vals[1, IPLabel- 1]-2*fac1*ff_vals[0, IPLabel- 1])

    return dFdt
