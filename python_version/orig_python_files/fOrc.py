from CalcBCmSW import CalcBCmSW
from CalcBCpSW import CalcBCpSW
from Calcf0 import Calcf0
import numpy as np

def fOrc(t, Start, TCoeffs, d, rmaxp1, Tot_Int_Pts, Gamma, Del_x, A, Shock_Flag, Exit_Pressure):
    #FORC evaluates ODE driver function f at l**{s}_{i}(t) at interior grd-pts
    #
    #   fOrc evaluates the ODE driver function at l**{s}_{i}(t) at the knot time
    #       t in subsubinterval j for each interior grid-point. ODE driver 
    #       function f is for 1D inviscid, compressible flow through a nozzle.
    #
    #   INPUTS: t = knot time value
    #           Start = starting time for subsubinterval j
    #           TCoeffs = d x rmaxp1 x Tot_Int_Pts array with Taylor 
    #                       Polynomial coefficients for l**{s}_{i}(t) in 
    #                       subsubinterval j at each interior grid-point
    #           d = number of components of f and l**{s}_{i}(t)
    #           rmaxp1 = number of terms/coefficients in a Taylor polynomial
    #           Tot_Int_Pts = number of interior grid-points
    #           Gamma = ratio of specific heats
    #           Del_x = distance between grid-points
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #           Shock_Flag = 0 (1) if shock wave absent (present)
    #           Exit_Pressure = pressure at nozzle exit
    #
    #   OUTPUT:
    #           f_Loc = d x Tot_Int_Pts array storing d components of ODE 
    #                   driver function f(U) at each interior grid-point at  
    #                   current value of primary flow variable U
    #
    #   Support functions: CalcBCmSW, CalcBCpSW, Calcf0

    # initialize parameter and array

    Tot_X_Pts = Tot_Int_Pts + 2 # number of grid-points
    delt = t - Start            # elapsed time from start of subsubinterval j
    lt = np.zeros((d,Tot_Int_Pts))   # stores l**{s}_{i}(t) at each interior 
                        #     grid-point

    Poly = np.zeros(rmaxp1)  # initialize to zero array storing Taylor 
                    #    polynomial coefficients for l(t(j-1), i) 
                    #    for given component and interior grid-point
                    
    U = np.zeros((d, Tot_X_Pts))  # array to store primary flow variables

    # evaluate l**{s}_{i}(t) at each interior grid-point, one component at time

    for ll in range(Tot_Int_Pts):    
        for m in range(d):
            # load Taylor polynomial coefficients for m-th component at
            #   interior grid-point ll in array Poly

            for column in range(rmaxp1):
                Poly[column] = TCoeffs[m, column, ll]

            # store value of m-th Taylor polynomial at elapsed time delt and 
            #       at interior grid-point ll

            lt[m, ll] = np.polyval(Poly, delt) 

    # assign U at each interior grid-point

    for ll in range(1, Tot_X_Pts - 1):
        IP_Label = ll - 1

        for m in range(d):
            U[m, ll] = lt[m, IP_Label]

    # assign U at boundary points using flow boundary conditions

    if (Shock_Flag == 0):
        U_Bvals = CalcBCmSW(U,A,Gamma,d,Tot_X_Pts)
    elif (Shock_Flag == 1):
        U_Bvals = CalcBCpSW(U,A,Gamma,d,Tot_X_Pts,Exit_Pressure)
    else:
        print('Unknown Shock_Flag value: ',Shock_Flag)

    for m in range(d):
        U[m, 0] = U_Bvals[m, 0]
        U[m, Tot_X_Pts - 1] = U_Bvals[m, 1]

    # evaluate f using Calcf0

    f_Loc = Calcf0(d, Tot_X_Pts, Tot_Int_Pts, Gamma, Del_x, U, A)

    return f_Loc