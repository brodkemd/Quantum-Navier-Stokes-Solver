
import numpy as np

def CalcBCpSW(U,A,Gamma,d,Tot_X_Pts,Exit_Pressure):
    #CALCBCPSW evaluates flow variables at nozzle boundaries-shock-wave present
    #   CalcBCpSW evaluates flow variables at boundaries of a nozzle. Flow is
    #       1D compressible, inviscid flow with shock-wave present.
    #   
    #   INPUTS:
    #           U = d x Tot_X_Pts array with current values of flow variables
    #                   at interior grid-points only. Values at boundary
    #                   grid-points are zero - will be updated in calling
    #                   program.
    #           A = 1 x Tot_X_Pts array with values of nozzle area at
    #                   grid-points.
    #           Gamma = ratio of specific heats
    #           d = number of flow variables at a grid-point
    #           Tot_X_Pts = number of grid-points
    #           Exit_Pressure = pressure at nozzle exit
    #
    #   OUTPUT:
    #           U_Bvals = d x 2 array with values of flow variable at nozzle
    #                       entrance and exit.

    U_Bvals = np.zeros((d,2))

    fac1 = 1/(Gamma -1)
    fac2 = Gamma/2

    GP_Nm1 = Tot_X_Pts -1
    GP_Nm2 = Tot_X_Pts -2

    # evaluate U at nozzle entrance:

    U_Bvals[0, 0] = A[0]
    U_Bvals[1, 0] = 2*U[1, 1] - U[1, 2]

    fac3 = (U_Bvals[1, 0]/U_Bvals[0, 0])**(2)

    U_Bvals[2, 0] = U_Bvals[0, 0]*(fac1 + fac2*fac3)

    # evaluate flow at nozzle exit:

    for k in range(d):
        U_Bvals[k, 1] = 2*U[k, GP_Nm1] - U[k, GP_Nm2]

    U_Bvals[2, 1] = Exit_Pressure*A[Tot_X_Pts- 1]*fac1 +( fac2*( U_Bvals[1, 1])**(2) )/U_Bvals[0, 1]

    return U_Bvals