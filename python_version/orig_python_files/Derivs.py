
from CalcFlux import CalcFlux
from CalcSource import CalcSource
from CalcFunc import CalcFunc
from Calc_d2Jdt2 import Calc_d2Jdt2
from Calc_d2Fdt2 import Calc_d2Fdt2
from Calc_dFdt import Calc_dFdt
from Calc_dffdt import Calc_dffdt
from Calc_dJdt import Calc_dJdt
from Calc_d2ffdt2 import Calc_d2ffdt2
from CalcfBvalsmSW import CalcfBvalsmSW
from CalcfBvalspSW import CalcfBvalspSW
from CalcdfdtBvalsmSW import CalcdfdtBvalsmSW
from CalcdfdtBvalspSW import CalcdfdtBvalspSW
from pprint import pprint
import numpy as np

def Derivs( d,r,InitVal,Del_x,Gamma,Tot_Int_Pts,Tot_X_Pts, A,Shock_Flag):
    #DERIVS evaluates ODE driver function f and its first r time derivatives
    #   Derivs evaluates the ODE driver function f(z) and its first r time
    #    derivatives at z = InitVal. 
    #
    #   INPUTS:
    #           d = number of components of ODE driver function
    #           r = maximum number of time derivatives to evaluate
    #           InitVal = d x Tot_X_Pts array storing flow variables at
    #                           all grid-points
    #           Del_x = separation between grid-points
    #           Gamma = ratio of specific heats
    #           Tot_Int_Pts = number of interior grid-points
    #                       = Tot_X_Pts - 2
    #           Tot_X_Pts = number of grid-points
    #           A = 1 x Tot_X_Pts array storing nozzle area at each grid-point
    #           Shock_Flag = 0(1) if shock-wave absent (present) in nozzle
    #
    #   OUTPUT: 
    #           ff = d x (r+1) x Tot_Int_Pts array where each page corresponds 
    #                   to a grid-point, and on each page the column indexes 
    #                   time derivatives, and row indexes the components of 
    #                   f(z). 
    #
    # Support functions: CalcFlux, CalcSource, CalcFunc, CalcfBvalsmSW,
    #                       CalcfBvalspSW, Calc_dFdt, Calc_dJdt, Calc_dffdt, 
    #                           CalcdfdtBvalsmSW, CalcdfdtBvalspSW,
    #                           Calc_d2Fdt2, Calc_d2Jdt2, Calc_d2ffdt2

    rmax = r+1     # max column label (store fcn and first r derivatives)

    # F = d x Tot_X_Pts array - stores fluxes contributing to ff at grid-points

    # J = 1 x Tot_Int_Pts array - stores source current contributing to ff
    #       at interior grid-points

    ff = np.zeros((d,rmax,Tot_Int_Pts)) # store i-th derivative values in column i
    #  for each component of driver function and
    #  interior grid-point


    # the following formulas were derived in my black research notebook entry
    #   for 03/08/2019.

    # Calculate ODE driver function ff at all interior grid-points
    #
    # 1.a Evaluate flow fluxes F which is a d x Tot_X_Pts array

    F = CalcFlux(InitVal,Gamma,d,Tot_X_Pts)


    # 1.b Evaluate source current J in momentum equation of motion which is
    #       a 1 x Tot_Int_Pts array

    J = CalcSource(InitVal,A,Gamma,Tot_X_Pts)


    # 1.c Now calculate ff(d,1, Tot_Int_Pts) which stores driver function 
    #       values. First calculate ff_vals = d x Tot_Int_Pts array which 
    #       stores ff values at all interior grid-point.

    ff_vals = CalcFunc(F,J,Del_x,d,Tot_Int_Pts)


    # 1.d Now store ff_vals in d x rmax x Tot_Int_Pts array ff

    for ll in range(Tot_Int_Pts):
        for k in range(d):
            #print("right =", ff_vals[k, ll])
            ff[k, 0, ll] = ff_vals[k, ll]


    # determine first (time) derivatives of ODE driver function at all
    #   interior grid-points

    #   2.a calculate "values" of ODE driver function at boundary points - I
    #        really mean value of dU(k,1l =1,Tot_X_Pts )/dt at boundary points
    #        which are found from the boundary conditions. I abuse notation
    #        and denote these derivatives by f(k,1) and f(k,Tot_X_Pts) which
    #        are stored in columns 1 and 2 of ff_Bvals.
    #
    #       ff_Bvals = d x 2 array storing "values" of ODE driver function
    #                       at nozzle entrance (exit) in column 1 (2).

    if Shock_Flag == 0:
        ff_Bvals = CalcfBvalsmSW(InitVal,Gamma,ff_vals,d,Tot_Int_Pts)

    elif Shock_Flag == 1:
        ff_Bvals = CalcfBvalspSW(InitVal,Gamma,ff_vals,d,Tot_Int_Pts)

    else:
        print('Unkown Shock_Flag value: ',Shock_Flag)

    #   2.b calculate first time derivative of flow fluxes at all grid-points
    #       
    #       dFdt = d x Tot_X_Pts array storing the flow flux first time
    #                               derivatives at all grid-points.

    dFdt = Calc_dFdt(InitVal,ff_vals,ff_Bvals,Gamma,d,Tot_X_Pts)

    #   2.c calculate first time derivative of flow source term at all
    #           interior grid-points
    #
    #       dJdt = 1 x Tot_Int_Pts array storing the flow source term
    #                               first time derivative at all interior
    #                               grid-points

    dJdt = Calc_dJdt(InitVal,ff_vals,A,Gamma,d,Tot_Int_Pts)


    #   2.d calculate dffdt_vals = d x Tot_Int_Pts array which 
    #       stores values of first time derivative of ff at all interior 
    #       grid-points.

    dffdt_vals = Calc_dffdt(dFdt, dJdt, Del_x, d, Tot_Int_Pts)


    #   2.e Now store dffdt_vals in column 2 of d x rmax x Tot_Int_Pts array 
    #       ff

    for ll in range(Tot_Int_Pts):
        for k in range(d):
            ff[k, 1, ll] = dffdt_vals[k, ll]



    # determine second (time) derivatives of ODE driver function at all
    #   interior grid-points

    #   3.a calculate "values" first time derivative of ODE driver function 
    #        at boundary points - I really mean value of 
    #        d**(2)U(k,1l =1,Tot_X_Pts )/dt**(2) at boundary points which are 
    #        found from the boundary conditions. I abuse notation and denote 
    #        these derivatives by df(k,1)/dt and df(k,Tot_X_Pts)/dt which are 
    #        stored in columns 1 and 2 of dffdt_Bvals.
    #
    #       dffdt_Bvals = d x 2 array storing "values" of first derivative of
    #                       ODE driver function at nozzle entrance (exit) in 
    #                       column 1 (2).

    if Shock_Flag == 0:
        dffdt_Bvals = CalcdfdtBvalsmSW(InitVal,Gamma,ff_Bvals,dffdt_vals,d,Tot_Int_Pts)

    elif Shock_Flag == 1:
        dffdt_Bvals = CalcdfdtBvalspSW(InitVal,Gamma,ff_Bvals,dffdt_vals,d, Tot_Int_Pts)

    else:
        print('Unknown Shock_Flag value: ', Shock_Flag)
    

    #   3.b calculate second time derivative of flow fluxes at all grid-points
    #       
    #       d2Fdt2 = d x Tot_X_Pts array storing the flow flux second time
    #                               derivatives at all grid-points.

    d2Fdt2 = Calc_d2Fdt2(InitVal,dffdt_vals,dffdt_Bvals,ff_vals,ff_Bvals, Gamma,d,Tot_X_Pts)


    #   3.c calculate second time derivative of flow source term at all
    #           interior grid-points
    #
    #       d2Jdt2 = 1 x Tot_Int_Pts array storing the flow source term
    #                               second time derivative at all interior
    #                               grid-points

    d2Jdt2 = Calc_d2Jdt2(InitVal,dffdt_vals,ff_vals,A,Gamma,d,Tot_Int_Pts)

    #   3.d calculate d2ffdt2_vals = d x Tot_Int_Pts array which 
    #       stores values of second time derivative of ff at all interior 
    #       grid-points.

    d2ffdt2_vals = Calc_d2ffdt2(d2Fdt2, d2Jdt2, Del_x, d, Tot_Int_Pts)


    #   3.e Now store d2ffdt2_vals in column 3 of d x rmax x Tot_Int_Pts array 
    #       ff
    
    '''
    print("f shape =", ff.shape)
    print("d shape =", d2ffdt2_vals.shape)
    print("d =", d)
    print("Tot_Int_Pts =", Tot_Int_Pts)
    '''
    for ll in range(Tot_Int_Pts):
        for k in range(d):
            ff[k, 2, ll] = d2ffdt2_vals[k, ll]

    #pprint(ff)
    #exit()

    return ff

