
import numpy as np
import sys


def Derivs(self, U):
        #DERIVS evaluates ODE driver function f and its first r time derivatives
      
        # F = d x Tot_X_Pts array - stores fluxes contributing to ff at grid-points
        # J = 1 x Tot_Int_Pts array - stores source current contributing to ff
        #       at interior grid-points
        #  for each component of driver function and
        #  interior grid-point
        ff = np.zeros((self.d,self.r+1,self.Tot_Int_Pts)) # store i-th derivative values in column i

        # the following formulas were derived in my black research notebook entry
        #   for 03/08/2019.

        # Calculate ODE driver function ff at all interior grid-points
        #
        # 1.a Evaluate flow fluxes F which is a d x Tot_X_Pts array
        F = self.CalcFlux(U)

        # 1.b Evaluate source current J in momentum equation of motion which is
        #       a 1 x Tot_Int_Pts array
        J = self.CalcSource(U)

        # 1.c Now calculate ff(d,1, Tot_Int_Pts) which stores driver function 
        #       values. First calculate ff_vals = d x Tot_Int_Pts array which 
        #       stores ff values at all interior grid-point.
        ff_vals = self.CalcFunc(F,J)

        # 1.d Now store ff_vals in d x rmax x Tot_Int_Pts array ff
        for ll in range(self.Tot_Int_Pts):
            for k in range(self.d):
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
        if self.Shock_Flag == 0: ff_Bvals = self.CalcfBvalsmSW(U, ff_vals)
        elif self.Shock_Flag == 1: ff_Bvals = self.CalcfBvalspSW(U, ff_vals)
        else: print('Unkown Shock_Flag value: ', self.Shock_Flag)

        #   2.b calculate first time derivative of flow fluxes at all grid-points
        #       
        #       dFdt = d x Tot_X_Pts array storing the flow flux first time
        #                               derivatives at all grid-points.
        dFdt = self.Calc_dFdt(U, ff_vals, ff_Bvals)
        
        #   2.c calculate first time derivative of flow source term at all
        #           interior grid-points
        #
        #       dJdt = 1 x Tot_Int_Pts array storing the flow source term
        #                               first time derivative at all interior
        #                               grid-points
        dJdt = self.Calc_dJdt(U, ff_vals)

        #   2.d calculate dffdt_vals = d x Tot_Int_Pts array which 
        #       stores values of first time derivative of ff at all interior 
        #       grid-points.
        dffdt_vals = self.Calc_dffdt(dFdt, dJdt)

        #   2.e Now store dffdt_vals in column 2 of d x rmax x Tot_Int_Pts array 
        #       ff
        for ll in range(self.Tot_Int_Pts):
            for k in range(self.d):
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
        if self.Shock_Flag == 0: dffdt_Bvals = self.CalcdfdtBvalsmSW(U, ff_Bvals, dffdt_vals)
        elif self.Shock_Flag == 1: dffdt_Bvals = self.CalcdfdtBvalspSW(U, ff_Bvals, dffdt_vals)
        else: print('Unknown Shock_Flag value: ', self.Shock_Flag)
        

        #   3.b calculate second time derivative of flow fluxes at all grid-points
        #       
        #       d2Fdt2 = d x Tot_X_Pts array storing the flow flux second time
        #                               derivatives at all grid-points.
        d2Fdt2 = self.Calc_d2Fdt2(U, dffdt_vals, dffdt_Bvals, ff_vals, ff_Bvals)


        #   3.c calculate second time derivative of flow source term at all
        #           interior grid-points
        #
        #       d2Jdt2 = 1 x Tot_Int_Pts array storing the flow source term
        #                               second time derivative at all interior
        #                               grid-points
        d2Jdt2 = self.Calc_d2Jdt2(U, dffdt_vals,ff_vals)

        #   3.d calculate d2ffdt2_vals = d x Tot_Int_Pts array which 
        #       stores values of second time derivative of ff at all interior 
        #       grid-points.
        d2ffdt2_vals = self.Calc_d2ffdt2(d2Fdt2, d2Jdt2)

        #   3.e Now store d2ffdt2_vals in column 3 of d x rmax x Tot_Int_Pts array 
        #       ff
        for ll in range(self.Tot_Int_Pts):
            for k in range(self.d):
                ff[k, 2, ll] = d2ffdt2_vals[k, ll]

        return ff


def BldTPoly(self, U):
        #BLDTPOLY calculates all Taylor polynomial coefficients for subint i
        
        # ll will store rmaxp1 Taylor poly. coeffs. for each component of ODE 
        # driver function, subsubinterval j, and interior grid-point.
        ll = np.zeros((int(self.d), int(self.r + 2), int(self.Tot_Int_Pts), self.N))   
        
        for j in range(self.N):
            # evaluate ODE driver func. and first rr time derivatives at InVal
            #    store in ff = dd x (rr+1) x Tot_Int_Pts array
            ff = self.Derivs(U)

            # for each subsubinterval j, calculate and store Taylor poly coeffs. 
            #  for each component of ODE driver function and interior grid-point  
            #  in array mat = dd x (rr+2) x Tot_Int_Pts 
            mat = self.BldTMat(ff)  

            ll[:, :, :, int(j)] = mat   # transfer Taylor polys coeffs for each component
                        # of driver function, subsubint j and interior
                        # grid-point

            
            U = self.NextInCond(U, mat) 

        # if change store values of ff at ithroat to send back to 
        #   main program on last iteration
        return ll, ff[:, :, int(self.ithroat) - 1]


def IntegrateGij(self, StoreLz, Start, i):
    #INTEGRATEGIJ integrates g_ij over subinterval i at each interior grd-pt
    # define array ti store integral result for subsubint j at each interior 
    #   grid-point
    IntegralValue = np.zeros((self.d, self.Tot_Int_Pts))   

    # define array to store integral result for entire subinterval i at each 
    #   interior grid-point
    Gint = np.zeros((self.d, self.Tot_Int_Pts))

    # Gij = d x Tot_Int_Pts x N array storing d components of g_ij at each 
    #   interior grid-point and N knot times for subsubinterval j
    # initialize parameters and arrays
    Gij = np.zeros((int(self.d), int(self.Tot_Int_Pts), self.N))

    #stor_time = []
    # loop over the subsubintervals j accumulate integral of driver function
    #   g_ij over subsubintervals.
    max_count = str(self.d * self.Tot_Int_Pts)
    for j in range(self.N):
        #start = time.time()
        print(f'Sub-interval: {j} / {self.N - 1}        ')

        # define array to store N knot times for sub-subinterval j
        if j == 0: t = np.linspace(Start, self.t[:, i][j], self.N)
        else: t = np.linspace(self.t[:, i][j-1], self.t[:, i][j], self.N)

        # evaluate f at N knot times for subsubinterval j and each interior
        # grid-point
        for k in range(self.N): Gij[:, :, k] = self.fOrc(t[k], t[0], StoreLz[:, :, :, j])

        # integrate Gij over subsubinterval j for each interior grid-point, 
        #  one component at a time, following basic approach in quantum 
        #  integration algorithm of Novak, J. Complexity vol. 17, 2-16
        #  (2001).
        #start_sub_for = time.time()
        count = 0
        for k in range(self.d):
            for ll in range(self.Tot_Int_Pts):
                count+=1
                print(f"{count} / {max_count}         ", end = "\r")
                GijVals = Gij[k, ll]
                # introduce gijVals which is a shifted and rescaled version 
                #   of GijVals that has values in range [0,1]
                # to that end, need to find max and min values of GijVals
                GijMax = np.max(GijVals)
                GijMin = np.min(GijVals)
                # will need the difference in these values
                DelGij = GijMax - GijMin
                
                # makes sure valid
                if DelGij <= self.Tolerance:
                    print('DelGij < Tolerance! Dividing by 0!')
                    exit(0)

                # use MeanOrc to evaluate mean of gijVals over N knot times
                #aTrue = self.MeanOrc(gijVals)
                omega = np.arcsin( np.sqrt( np.mean((GijVals - GijMin)/DelGij) ) ) / np.pi

                # use QAmpEst to estimate mean of gijVals over subsubint j
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
                # need to undo shift and rescaling to get integral of GijVals
                #   formula used is derived in Supplementary Material for 
                #   my paper describing this work.
                #temp_start = time.time()
                IntegralValue[k, ll] = self.hbar * (self.QAmpEst(omega) * DelGij + GijMin)
                #stor_time.append(time.time() - temp_start)
            

        #print("average time spent on QAmpEst=", np.mean(stor_time))
        #print("total time spent on QAmpEst=", self.d * self.Tot_Int_Pts * np.mean(stor_time))
        #print("time in sub for loop=", time.time() - start_sub_for)

        # add IntegralValue for subsubinterval j to integral over previous
        #      subsubintervals - at loop's end contains integral of driver
        #      function over subinterval i (to within factor of 
        #      1/N - see below)
        Gint = Gint + IntegralValue

        #dur = time.time() - start
        #print("time =", dur)
        #print("predicted time =", 16 * 256 * dur / 3600, "hours")
        #exit()
        sys.stdout.write("\033[F") 

    # Eq. (28) in Kacewicz requires Gint above to be divided by N ( = ml)
    #  other division by N included in calculation of mean by MeanOrc
    return Gint/self.N


def IntegrateODE(self):
    #INTEGRATEODE numerically integrates ODE for 1D Navier-Stokes flow
    print("Running...\nPress ctrl-C to stop the execution\n")
    
    # Begin loop over the subintervals i result is approximate solution z(t).
    for i in range(int(self.n)):
        print('Interval:', i)
        #build Taylor polynomials l**{s}_{i}(t)for subinterval i at all
        #   interior grid-points store polynomial coefficients in 
        #       StoreLz(d,r+2,Tot_Int_Pts,N)
        #   ff_throat is a d x (r+1) array storing the values of residual and
        #       its first r time derivatives at the nozzle throat at the end 
        #       of subinterval i
        StoreLz, ff_throat = self.BldTPoly(self.InitVal)

        # store values of ff_throat for subinterval i for easy of writing to
        #   files
        self.ff0_throat[:, i] = np.abs(ff_throat[:, 0])
        self.ff1_throat[:, i] = np.abs(ff_throat[:, 1])
        self.ff2_throat[:, i] = np.abs(ff_throat[:, 2])

        # define Start time for subinterval i
        if i == 0: Start = self.a
        else: Start = self.t[self.N - 1, i - 1]

        #  gInt = d x Tot_Int_Pts array storing integral of each component of
        #               g_ij over subinterval i for each interior grid-point
        gInt = self.IntegrateGij(StoreLz, Start, i)

        # add gInt for subint i to InitVal to get InitVal for 
        #   subinterval i + 1 at each interior grid-point
        for ll in range(1, (self.Tot_Int_Pts+1)):
            for m in range(self.d): self.InitVal[m, ll] = self.InitVal[m, ll] + gInt[m, ll - 1]

        # use flow boundary conditions to determine InitVal for
        #   subinterval i + 1 at two boundary points 1 and Tot_X_Pts. Store
        #   them in InitVal_BVals
        if self.Shock_Flag == 0: InitVal_BVals = self.CalcBCmSW(self.InitVal)
        elif self.Shock_Flag == 1: InitVal_BVals = self.CalcBCpSW(self.InitVal)
        else: print('Unknown Shock_Flag value: ', self.Shock_Flag)

        # transfer boundary values from InitVal_Bvals to InitVal. Resulting 
        #       InitVal then contains initial values of primary flow 
        #       variable U for subinterval i + 1 at each grid-point.
        for m in range(self.d):
            self.InitVal[m, 0] = InitVal_BVals[m, 0]
            self.InitVal[m, self.Tot_X_Pts - 1] = InitVal_BVals[m, 1]

        # store U2 values at start of subinterval i + 1 at all gridpoints
        for gridpt in range(self.Tot_X_Pts): self.U2[gridpt, i+1] = self.InitVal[1, gridpt]

        # evaluate physical flow variables from InitVal for next subinterval
        self.Calc_FlowVarResults()

        # write computational results to files for eventual plotting:
        self.WriteResults()
        
        # moving up a line
        sys.stdout.write("\033[F")

    print("\nDone...\n")