from fOrc import fOrc
import numpy as np

def FuncOrc(t, TCoeffs, d, rmaxp1, N, rho, Tot_Int_Pts, A, Gamma, Del_x, Shock_Flag, Exit_Pressure):
    #FUNCORC evaluates g_ij[u] at N knot times in subsubinterval j
    #
    #   FuncOrc evaluates g_ij[u] at N knot times for subsubinterval j
    #
    #   INPUTS: t = array of N knot time values in subsubinterval j
    #           TCoeffs = d x rmaxp1 x Tot_Int_Pts array of Taylor polynomial 
    #                       coefficients for d components of l**{s}_{i}(t) in 
    #                       subsubinterval j at each interior grid-point
    #           d = number of components of g_ij[u]
    #           rmaxp1 = number of terms/coefficients in a Taylor polynomial
    #           N = number of knot times in subsubinterval j
    #           rho = basic Holder class parameter.
    #           Tot_Int_Pts = number of interior grid-points
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #           Gamma = ratio of specific heats
    #           Del_x = distance between grid-points
    #           Shock_Flag = 0 (1) if shock wave absent (present)
    #           Exit_Pressure = pressure at nozzle exit
    #
    #   OUTPUT:
    #           Gij = d x Tot_Int_Pts x N array storing value of d components 
    #                       of g_ij at each interior grid-point and at the N
    #                       knot times for subsubinterval j
    #
    #   Support function: fOrc

    # initialize parameters and arrays

    r = rmaxp1 - 2     # number of derivatives of f allowed
    q = r + rho        # exponent of hbar in Gij formula 

    # f stores d components of ODE driver function f at each interior 
    #    grid-point and N knot times for subsubinterval j

    f = np.zeros((int(d), int(Tot_Int_Pts), int(N)))   

    # evaluate f at N knot times for subsubinterval j and each interior
    # grid-point

    for k in range(int(N)):
        f[:, :, k] = fOrc(t[k], t[0], TCoeffs, d, rmaxp1,Tot_Int_Pts,Gamma,Del_x,A, Shock_Flag, Exit_Pressure)


    # assign values of Gij at N knot times t for subsubinterval j and each 
    #       interior grid-point

    Gij = f

    return Gij