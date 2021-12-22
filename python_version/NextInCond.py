from CalcBCmSW import CalcBCmSW
from CalcBCpSW import CalcBCpSW
import numpy as np

def NextInCond( mmat, InitVal, hbar, d, r, A, Gamma, Tot_Int_Pts, Tot_X_Pts, Shock_Flag, Exit_Pressure):
    #NEXTINCOND determines initial condition for next subsubinterval 
    #   NextInCond evaluates the Taylor polynomials for a given subsubinterval
    #   at the subsubinterval's greatest value and uses this result as the
    #   initial value for the next subsubinterval.
    #
    #   INPUTS: mmat = d x (r+2) x Tot_Int_Pts array of Taylor polynomial 
    #                   coefficients for each driver function component and
    #                   interior grid-point
    #           InitVal = d x Tot_X_Pts array with initial values of primary
    #                       flow variables for current subsubinterval
    #           hbar = width of subsubinterval
    #           d = number of components of InitVal
    #           r + 1 = degree of Taylor polynomials
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #           Gamma = ratio of specific heats
    #           Tot_Int_Pts = number of interior grid-points
    #           Tot_X_Pts = number of grid-points
    #           Shock_Flag = 0 (1) if shock wave present (absent) in divergent
    #                           part of nozzle
    #           Exit_Pressure = pressure at nozzle exit
    #
    #   OUTPUT:
    #           NextInVal = d x Tot_X_Pts array with initial condition values 
    #                           for primary flow variables for next
    #                           subsubinterval
    #
    # Support functions: CalcBCmSW, CalcBCpSW

    rmax = r + 1       # degree of Taylor Polynomials
    rmaxp1 = rmax + 1  # number of terms in Taylor polynomials

    TylrPoly = np.zeros(rmaxp1) # array to store coefficients for a Taylor
        #   polynomial
        
    NextInVal = np.zeros((d,Tot_X_Pts))     # array to store initial conditions 
                #   for next subsubinterval 

    # calculate initial condition for next subsubinterval for the
    #       interior grid-points    ( 2 <= ll <= TotXPtm1 )


    for ll in range(Tot_Int_Pts):
        GPLabel = ll + 1

        for pp in range(d):
            for mm in range(rmax):
                TylrPoly[mm] = mmat[pp, mm, ll]

            TylrPoly[rmaxp1- 1] = InitVal[pp, GPLabel]

            NextInVal[pp, GPLabel] = np.polyval(TylrPoly,hbar)


    # use flow boundary conditions to determine new initial conditions for
    #       boundary grid-points

    if (Shock_Flag == 0):
        U_Bvals = CalcBCmSW(InitVal,A,Gamma,d,Tot_X_Pts)
    elif (Shock_Flag == 1):
        U_Bvals = CalcBCpSW(InitVal,A,Gamma,d,Tot_X_Pts,Exit_Pressure)
    else:
        print(' Unknown Shock_Flag value: ', Shock_Flag)


    # store boundary values in NextInVal

    for p in range(d):
        NextInVal[p, 0] = U_Bvals[p, 0]      # values at nozzle entrance
        NextInVal[p, Tot_X_Pts- 1] = U_Bvals[p, 1]  # values at nozzle exit

    return NextInVal