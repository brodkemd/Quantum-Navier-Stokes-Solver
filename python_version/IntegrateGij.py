from FuncOrc import FuncOrc
from MeanOrc import MeanOrc
from QAmpEst import QAmpEst
import numpy as np

def IntegrateGij(StoreLz,StoreTimes4i,Start, d,r,N,delta1,hbar,rho, Tot_Int_Pts,A,Gamma, Del_x,Shock_Flag,i,Exit_Pressure):
    #INTEGRATEGIJ integrates g_ij over subinterval i at each interior grd-pt
    #   
    #   IntegrateGij integrates g_ij(l**{s}_{i}(t)) over subinterval i at each
    #       interior grid-point. g_ij is defined in Kacewicz, J. Complexity 
    #       vol. 22, 676-690 (2006). HERE g_ij is the ODE driver function f!!!
    #
    #   INPUTS: StoreLz = d x rmaxp1 x Tot_Int_Pts x N array with Taylor 
    #                       polynomial coefficients for l**{s}_{i}(t) for each
    #                       interior grid-point, component, & subsubinterval j
    #           StoreTimes4i = 1 x N array with times defining subsubintervals
    #                           in subinterval i. For j >= 2, start (end) time
    #                           is t(j-1,i) (t(j,i)) for i >= 2, j = 1, start
    #                           time is t(N, i-1) and for i = 1, j = 1, start
    #                           time is a which is the initial time for the ODE
    #                           solution
    #           Start = start time of subinterval i
    #           d = number of components of g_ij
    #           r + 2 = rmaxp1 = number of terms/coefficients in Taylor
    #                               polynomial l**{s}_{i}(t)
    #           N = number of subsubintervals and (as required in Kacewicz
    #                quantum algorithm) number of knot times in a
    #                subsubinterval
    #           delta1 = probability quantum integration result violates its
    #                       error upper bound
    #           hbar = width of each subsubinterval
    #           rho = basic Holder class parameter.
    #           Tot_Int_Pts = number of interior grid-points
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #           Gamma = ratio of specific heats
    #           Del_x = distance between grid-points
    #           Shock_Flag = 0 (1) if shock wave absent (present)
    #           Exit_Pressure = pressure at nozzle exit
    #
    #   OUTPUT:
    #           Gint = d x Tot_Int_Pts array storing d components of integral
    #                       of Gij over subinterval i for each interior
    #                       grid-point
    #
    #   Support functions: FuncOrc, MeanOrc, QAmpEst

    # initialize parameter and arrays        print("InVal =", InVal)

    rmaxp1 = r + 2

    Tolerance = 10**(-1*12)  # used to test for division by zero below

    # define array ti store integral result for subsubint j at each interior 
    #   grid-point

    IntegralValue = np.zeros((d, Tot_Int_Pts))   

    # define array to store integral result for entire subinterval i at each 
    #   interior grid-point

    Gint = np.zeros((d, Tot_Int_Pts))

    # loop over the subsubintervals j accumulate integral of driver function
    #   g_ij over subsubintervals.

    for j in range(int(N)):
        print('In subinterval i =', i, '   starting subsubinterval j =',j)

        # define array to store N knot times for sub-subinterval j
        
        #print("start =", Start)
        #print("StoreTimes4i[j] =", StoreTimes4i[j])

        if j == 0:
            t = np.linspace(Start, StoreTimes4i[j], int(N))
        else:
            t = np.linspace(StoreTimes4i[j-1], StoreTimes4i[j], int(N))

        #print("t =", t)
        #print(t.shape)


        # Gij = d x Tot_Int_Pts x N array storing d components of g_ij at each 
        #   interior grid-point and N knot times for subsubinterval j


        Gij = FuncOrc(t,StoreLz[:, :, :, j],d,rmaxp1,N, rho, Tot_Int_Pts, A, Gamma, Del_x, Shock_Flag, Exit_Pressure)

        #print("Gij =", Gij)

        # GijVals stores values of Gij (viz. driver function f) at N knot times
        #   for sub-subinterval j for a given component k of Gij and interior
        #   grid-point ll

        GijVals = np.zeros((int(N)))

        # integrate Gij over subsubinterval j for each interior grid-point, 
        #  one component at a time, following basic approach in quantum 
        #  integration algorithm of Novak, J. Complexity vol. 17, 2-16
        #  (2001).

        for k in range(d):
            for ll in range(Tot_Int_Pts):
                for knot in range(int(N)):
                    GijVals[knot] = Gij[k, ll, knot]

                # introduce gijVals which is a shifted and rescaled version 
                #   of GijVals that has values in range [0,1]

                # to that end, need to find max and min values of GijVals

                GijMax = np.max(GijVals)
                GijMin = np.min(GijVals)

                # will need the difference in these values

                DelGij = GijMax - GijMin

                # test whether DelGij is small number

                if DelGij > Tolerance:
                    # now define gijVals

                    gijVals = (GijVals - GijMin)/DelGij

                    # use MeanOrc to evaluate mean of gijVals over N knot times

                    aTrue = MeanOrc(gijVals, t, N)

                    omega = np.asin(np.sqrt(aTrue))/np.pi

                    # use QAmpEst to estimate mean of gijVals over subsubint j

                    aEstimate = QAmpEst(N, delta1, omega)

                    # NOTE: aEstimate must be multiplied by hbar for mean 
                    #           estimate to approximate integral (this 
                    #           restores the delta-t)

                    # add contributions to integral of gijVals for 
                    #   subsubinterval j.
                    #
                    #   NOTE: (i) this integral is for gijVals, not GijVals - 
                    #               will correct shortly
                    #         (ii) aEstimate must be multiplied by hbar as 
                    #               noted above

                    IntegralValue[k,ll] = hbar*aEstimate

                elif DelGij <= Tolerance:
                    IntegralValue[k, ll] = 0.0

                    print('DelGij < Tolerance! Beware dividing by 0!')
                    print("GijMin =", GijMin)
                    print("GijMax =",GijMax)
                    print("DelGij =",DelGij)
                    input('   Press any key to continue calculation...')


                # need to undo shift and rescaling to get integral of GijVals
                #   formula used is derived in Supplementary Material for 
                #   my paper describing this work.

                IntegralValue[k,ll] = IntegralValue[k,ll]*DelGij + hbar*GijMin

        # add IntegralValue for subsubinterval j to integral over previous
        #      subsubintervals - at loop's end contains integral of driver
        #      function over subinterval i (to within factor of 
        #      1/N - see below)

        Gint = Gint + IntegralValue


    # Eq. (28) in Kacewicz requires Gint above to be divided by N ( = ml)
    #  other division by N included in calculation of mean by MeanOrc

    Gint = Gint/N

    return Gint
