
import numpy as np

def CalcdfdtBvalspSW(U,Gamma,ff_Bvals,dffdt_vals,d, Tot_Int_Pts):
    #CALCDFDTBVALSPSW assigns dff/dt boundary values - shock-wave present
    #   CalcdfdtBvalspSW uses flow boundary conditions to assign values to
    #       first time derivative of ODE driver function at nozzle entrance
    #       and exit when shock-wave present. Formulas below are based on 
    #       analysis in 7-23-2019 entry of my research notebook.
    #
    #   INPUTS:
    #           U = d x Tot_X_Pts array storing primary flow variables at all
    #                   grid-points.
    #           Gamma = ratio of specific heats
    #           ff_Bvals = d x 2 array storing "values" of ODE driver function
    #                           at nozzle entrance and exit.
    #           dffdt_vals = d x Tot_Int_Pts array storing values of first 
    #                           time derivative of ODE driver function at all
    #                           interior grid-points.
    #           d = number of primary flow variables at each grid-point
    #           Tot_Int_Pts = number of interior grid-points
    #
    #   OUTPUT:
    #           dffdt_Bvals = d x 2 array storing "values" of first time
    #                               derivative of ODE driver function at
    #                               nozzle entrance and exit.

    dffdt_Bvals = np.zeros((d,2))

    Fac1 = U[1, 0]/U[0, 0]
    Fac2 = (ff_Bvals[1, 0])**(2)/U[0, 0]

    # nozzle entrance calculations

    dffdt_Bvals[0, 0] = 0

    dffdt_Bvals[1, 0] = 2*dffdt_vals[1, 0] - dffdt_vals[1, 1]

    dffdt_Bvals[2, 0] = Gamma*( (Fac1*dffdt_Bvals[1, 0]) + Fac2 )

    # nozzle exit calculations

    Tot_X_Pts = Tot_Int_Pts + 2

    Fac0 = Gamma/2

    Fac3 = U[1, Tot_X_Pts- 1]/U[0, Tot_X_Pts - 1]

    Fac4 = Fac3**(2)
    Fac5 = 1/U[0, Tot_X_Pts- 1]
    Fac6 = Fac3*Fac5
    Fac7 = Fac4*Fac5

    for k in range(d):
        dffdt_Bvals[k, 1] = 2*dffdt_vals[k, Tot_Int_Pts- 1]- dffdt_vals[k, (Tot_Int_Pts - 2)]


    dffdt_Bvals[2, 1] = Fac0*( Fac3*( 2*dffdt_Bvals[1, 1]- Fac3*dffdt_Bvals[0, 1] ) +2*Fac7*(ff_Bvals[0, 1])**(2) +2*Fac5*(ff_Bvals[1, 1])**(2)-4*Fac6*ff_Bvals[0, 1]*ff_Bvals[1, 1] )

    return dffdt_Bvals
