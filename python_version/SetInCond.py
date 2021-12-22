import numpy as np
import random

def SetInCond(Shock_Flag, In_Mass_Flow_E,Gamma,x, Del_x,A,d,Mrho_E,Temp_E, ICMFlowErrScale, ICrhoErrScale,ICtempErrScale):
    #SETINCOND assigns initial condition for flow variables U
    #   
    #   INPUT: 
    #           Shock_Flag = 0 (1) if shock wave absent (present)
    #           In_Mass_Flow_E = mass flow value for exact solution of 
    #                               subsonic to supersonic flow w/o shock-wave 
    #           Gamma = ratio of specific heats C_p/C_v
    #           x = 1 x Tot_Pts array with x grid-point locations
    #           Del_x = spatial separation between grid-points
    #           A = 1 x Tot_Pts array with nozzle area at grid-points
    #           d = number of flow variables at each grid-point
    #           Mrho_E = 1 x Tot_X_Pts array storing exact values of mass
    #                       density at each grid-point
    #           Temp_E = 1 x Tot_X_Pts array storing exact values of 
    #                       temperature at each grid-point
    #           ICMFlowErrScale = scale of random shift to initial mass flow 
    #                               away from exact results 
    #           ICrhoErrScale = scale of random shift to initial mass density
    #                               away from exact results (no shift at nozzle  
    #                               throat or inlet)
    #           ICtempErrScale = scale of random shift to initial temperature
    #                               away from exact results (no shift at nozzle  
    #                               throat or inlet) 

    #
    #   OUTPUT:
    #           U = d x Tot_Pts array of flow variables
    #           Del_t = time-step size (determined by CFL condition)
    #           In_Mass_Flow = initial mass flow rate used in simulation:
    #                               includes random shift

    # determine number of grid-points

    Tot_Pts = len(x)

    # determine coefficients and set Tiny to small value for later use

    Coef1 = 1 /(Gamma -1)
    Coef2 = Gamma/2

    ithroat = (Tot_Pts + 1)/2

    # initialize arrays

    Mrho = np.zeros((Tot_Pts))    # will store shifted mass density
    Temp = np.zeros((Tot_Pts))    # will store shifted temperature
    U = np.zeros((d,Tot_Pts))       # will store flow variables
    V = np.zeros((Tot_Pts))       # will store velocity
    Loc_TSteps = np.zeros((Tot_Pts - 2))    # will store local time-step

    # add random shift to In_Mass_Flow_E

    In_Mass_Flow = In_Mass_Flow_E*(1 + (-1 * ICMFlowErrScale +2*ICMFlowErrScale*random.random()))

    # begin loop over grid-points - initialize mass density and temperature

    # Introduce a small shift

    #Tiny = 10**(-8)
    #print("Shock flag =", Shock_Flag)
    for i in range(int(Tot_Pts)):
        #print("I =", i)
        if Shock_Flag == 0:
            if ( (i == 0) or (i == ithroat - 1) ):
                Mrho[i] = Mrho_E[i]
                Temp[i] = Temp_E[i]
            elif ((i != 0) and (i != ithroat - 1) ):
                Mrho[i] = Mrho_E[i]*( 1 + (-1 *ICrhoErrScale + 2*ICrhoErrScale*random.random()) )

                if Mrho[i] > 1:
                    Mrho[i] = 1

                Temp[i] = Temp_E[i]*( 1 + (-1 * ICtempErrScale + 2*ICtempErrScale*random.random()) )

                if Temp[i] > 1:
                    Temp[i] = 1

        elif Shock_Flag == 1:
            if ((i == 0) or (i == ithroat - 1) ):
                Mrho[i] = Mrho_E[i]
                Temp[i] = Temp_E[i]
            elif ( (i != 0) and (i != ithroat - 1) ):
                Mrho[i] = Mrho_E[i]*( 1 + (-1 * ICrhoErrScale + 2*ICrhoErrScale*random.random()) )

                if Mrho[i] > 1:
                    Mrho[i] = 1

                Temp[i] = Temp_E[i]*( 1 + (-1 * ICtempErrScale + 2*ICtempErrScale*random.random()) )

                if Temp[i] > 1:
                    Temp[i] = 1
        # assign initial condition to flow variables

        U[0, i] = Mrho[i]*A[i]
        U[1, i] = In_Mass_Flow
        U[2, i] = U[0, i]*( Coef1*Temp[i] + Coef2*(U[1, i]/U[0, i])**(2) )

        #print("U[0, i] =", U[0, i])
        #print("U[1, i] =", U[1, i])
        #print("U[2, i] =", U[2, i])
        


    # determine time-step Del_t using Courant-Friedrichs-Levy (CFL) stability
    #   condition with C = 0.5

    C = 0.5

    for i in range(1, (Tot_Pts - 1)):
        V[i] = U[1, i]/U[0, i]
        Loc_TSteps[i - 2] = (C*Del_x)/(np.sqrt(Temp[i]) + V[i])

    # Del_t is the min value in Loc_TSteps

    Delta_t = np.min(Loc_TSteps)

    return U, Delta_t,In_Mass_Flow

