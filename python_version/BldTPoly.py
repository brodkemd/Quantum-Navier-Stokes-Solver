from Derivs import Derivs 
from BldTMat import BldTMat
from NextInCond import NextInCond
import numpy as np

def BldTPoly( dd,nn,NN,hb,rr,InVal,Del_x, Gamma,Tot_Int_Pts, Tot_X_Pts, A, Shock_Flag, Exit_Pressure, ithroat):
    #BLDTPOLY calculates all Taylor polynomial coefficients for subint i
    #   BldTPoly calculates the Taylor polynomial coefficients for each
    #   subsubinterval j in the subinterval i.
    #   INPUTS: dd = number of components of ODE solution
    #           nn = number of subintervals
    #           NN = number of subsubintervals j in a given subinterval i
    #           hb = width of a subsubinterval
    #           rr + 1 = degree of Taylor polynomials
    #           InVal = d x Tot_X_Pts array with initial condition for Taylor 
    #                       polynomials in subint i
    #           Del_x = separation between grid-points
    #           Gamma = ratio of specific heats
    #           Tot_Int_Pts = number of interior grid-points
    #           Tot_X_Pts = number of grid-points
    #           A = 1 x Tot_X_Pts array with cross-sectional nozzle area at
    #                   each grid-point
    #           Shock_Flag = 0 (1) if shock wave present (absent) in divergent
    #                           part of nozzle
    #           Exit_Pressure = pressure at nozzle exit
    #           ithroat = grid-point index for nozzle throat
    #
    #   OUTPUT:
    #           ll = dd x (rr + 2) x Tot_Int_Pts x NN array storing the
    #                       Taylor polynomial coefficients for each component
    #                       of ODE driver function at each interior grid-point
    #                       and subsubinterval
    #
    #           ff_throat = dd x (rr + 1) array that stores values of ODE
    #                           driver function and its first rr time
    #                           derivatives at the last subsubinterval NN
    #
    # Support functions: Derivs, BldTMat, NextInCond

    rmaxp1 = rr + 2    # number of terms in Taylor polynomials
    #print("rmaxp1 =", rmaxp1)
    # ll will store rmaxp1 Taylor poly. coeffs. for each component of ODE 
    # driver function, subsubinterval j, and interior grid-point.

    ll = np.zeros((int(dd),int(rmaxp1),int(Tot_Int_Pts),int(NN)))   
    #print(ll.shape)

    for j in range(int(NN)):
        # evaluate ODE driver func. and first rr time derivatives at InVal
        #    store in ff = dd x (rr+1) x Tot_Int_Pts array

        ff = Derivs(dd,rr,InVal,Del_x,Gamma,Tot_Int_Pts,Tot_X_Pts,A, Shock_Flag)

        # if j == NN store values of ff at ithroat to send back to 
        #   main program

        if j == NN - 1:
            ff_throat = ff[:, :, int(ithroat)- 1]

        # for each subsubinterval j, calculate and store Taylor poly coeffs. 
        #  for each component of ODE driver function and interior grid-point  
        #  in array mat = dd x (rr+2) x Tot_Int_Pts 

        mat = BldTMat(ff,dd,rr,InVal,Tot_Int_Pts,Tot_X_Pts)  

        ll[:, :, :, int(j)] = mat   # transfer Taylor polys coeffs for each component
                    # of driver function, subsubint j and interior
                    # grid-point

        
        #print(ll.shape)
        #exit()
        InVal = NextInCond(mat,InVal,hb,dd,rr,A,Gamma,Tot_Int_Pts, Tot_X_Pts,Shock_Flag,Exit_Pressure) 

    return ll, ff_throat

