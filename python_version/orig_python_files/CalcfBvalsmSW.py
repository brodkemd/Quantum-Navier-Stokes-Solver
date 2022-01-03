
import numpy as np

def CalcfBvalsmSW(U,Gamma,ff_vals,d,Tot_Int_Pts):
    #CALCFBVALSMSW assigns ODE driver function boundary values - no shock-wave
    #   CalcfBvalsmSW uses flow boundary conditions to assign "values" to ODE
    #       driver function at nozzle entrance/exit when no shock-wave is 
    #       present in nozzle. Actually calculating value of 
    #       dU(k,ll = 1, Tot_X_Pts)/dt at boundary points using boundary 
    #       conditions. Formulas below based on analysis in 4-23-2019 entry 
    #       in my research notebook.
    #
    #   INPUTS:
    #           U = d x Tot_X_Pts array storing flow variables at grid-points
    #           Gamma = ratio of specific heats
    #           ff_vals = d x Tot_Int_Pts array storing ODE driver function
    #                           values at interior grid-points
    #           d = number of components of ODE driver function
    #           Tot_Int_Pts = number of interior grid-points
    #
    #   OUTPUT:
    #           ff_Bvals = d x 2 array storing ODE driver function values at
    #                           nozzle entrance and exit

    ff_Bvals = np.zeros((d,2))

    Fac3 = U[1, 0]/U[0, 0]

    # calculate ff_Bvals at nozzle entrance

    ff_Bvals[0, 0] = 0
    ff_Bvals[1, 0] = 2*ff_vals[1, 0] - ff_vals[1, 1]
    ff_Bvals[2, 0] = Gamma*Fac3*ff_Bvals[1, 0]

    # calculate use boundary conditions on flow variables to determine
    #   ff_Bvals at nozzle exit

    for k in range(3):
        ff_Bvals[k, 1] = 2*ff_vals[k, Tot_Int_Pts] - ff_vals[k, (Tot_Int_Pts-1)]


    return ff_Bvals