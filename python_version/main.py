import time, warnings, random
import numpy as np

"""
how to see all varaibles in matlab
myvals = whos;
for n = 1:length(myvals)
      myvals(n).name
      eval(myvals(n).name)
end

in python 
dir()

"""

warnings.filterwarnings("ignore")      

class ns_q:
    a = 0
    Tot_TSteps = 1400
    in_n = 16
    d = 3
    r = 2
    err1 = 0.005
    delta = 0.005
    rho = 1
    x_min = 0
    x_max =3.0
    Tot_X_Pts = 61
    Gamma = 1.4
    Exit_Pressure = 0.6784
    Shock_Flag = 1

    # In_Mass_Flow set to exact value for sub-to-supersonic flow without
    #   shock-wave
    In_Mass_Flow = 0.579

    # ICMFlowErrScale less than or about 0.01 which is size of shift used by
    #   Anderson
    ICMFlowErrScale = 0.00

    # ICrhoErrScale is 10# (2#) of minimum density in steady-state solution
    #   when shockwave absent (present)
    ICrhoErrScale = 0.02

    # ICtempErrScale is 2# (1#) of minimum temperature in steady-state 
    #   solution when shockwave absent (present)
    ICtempErrScale = 0.01      


    def check(self, items_in_list):
        check_l = []
        for item in dir(self):
            if "__" not in item and not callable(eval("self." +item)):
                check_l.append(item)
        
        items_in_list = items_in_list.split(", ")

        for item in items_in_list:
            if item not in check_l:
                raise ValueError(item + " not in variables")


    def show(self, val=True):
        for item in sorted(dir(self), key=str.lower):
            if "__" not in item and not callable(eval("self." +item)):
                if val: print(item, "=", eval("self." + item))
                else: print(item)


    def __init__(self):
        start = time.time()

        self.InitCalcParms()

        # integrate ODE
        #self.IntegrateODE()

        '''
        # stop the timer for the calculation

        runtime = (start - time.time())/60
        timepersubint = runtime/n

        # write computational results to files for eventual plotting:

        #WriteResults(n, Tot_X_Pts, d, U2, Mach_D, Mach_E, Mrho_D, Mrho_E, Press_D, Press_E, Temp_D, Temp_E, Rel_MachErr, Rel_MrhoErr, Rel_PressErr, Rel_TempErr, AvRelTempErr, AvPlusSDevRelTempErr, AvMinusSDevRelTempErr, AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr, AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr, AvRelPressErr, AvPlusSDevRelPressErr, AvMinusSDevRelPressErr, AvU2, ff0_throat, ff1_throat, ff2_throat)
        
        message = 'QNavStokes_solvr has finished results written to files.'

        print(message)
        '''
        print('Program runtime (minutes) = ', (start - time.time())/60)
        
        #timepersubint = runtime/n
        #print('Program runtime per subinterval(minutes) = ', timepersubint)

        self.show(val = False)


    def InitCalcParms(self):
        # set up x-grid for problem
        self.x = np.linspace(self.x_min, self.x_max, self.Tot_X_Pts)
        self.Del_x = self.x[1] - self.x[0]

        # calculate nozzle area at each grid-point
        self.Calc_Noz_Area()

        # calculate exact steady-state solution of primary flow variables: Mach
        #   number, mass density, pressure, and temperature
        if self.Shock_Flag == 0: self.Calc_ExactResultsmSW()
        elif self.Shock_Flag == 1: self.Calc_ExactResultspSW()
        else: print('Unknown Shock_Flag value: ', self.Shock_Flag)

        # Initialize initial condition InitVal of ODE solution. Delta_t is maximum
        #   time-step satisfying Courant-Friedrichs-Levy (CFL) numerical stability
        #   condition. InitVal is d x Tot_X_Pts array storing initial 
        #   computational flow variables U(d,Tot_X_Pts)
        self.SetInCond()

        # calculate final time b
        b = self.Tot_TSteps*self.Delta_t

        # initialize parameters for quantum ODE algorithm
        self.InitParms()
        self.n = self.final_n

        # Partition time interval [a, b] into n^(k) sub-subintervals by introducing
        #   N - 1 intermediate times t(j,i), where t(j,i) is largest time in
        #   sub-subinterval j in subinterval i. The smallest time in
        #   sub-subinterval j is t(j-1,i). Note that the initial time t = a is NOT
        #   stored in N x n array t
        self.IPrtn(b) # b is a local so it is just passed in

        # Taylor polynomials only needed at interior grid-points so define
        self.Tot_Int_Pts = self.Tot_X_Pts - 2

        # initialize array U2 to store calculated mass flow rate at beginning of
        #   each subinterval i and also at time at end of ODE integration.
        #       U2 = Tot_X_Pts x (n+1) array

        self.U2_in = np.zeros((int(self.Tot_X_Pts), int(self.n+1)))

        # use initial condition to assign U2 at start of subinterval i = 1

        for gridpt in range(int(self.Tot_X_Pts)):
            self.U2_in[gridpt, 0] = self.InitVal[1, gridpt]

        # define ithroat as grid-point index of nozzle throat

        self.ithroat = (self.Tot_X_Pts + 1)/2

        # initialize arrays to store value of driver function f(z(t)) and its first
        #   r time derivatives at nozzle throat at end of each subinterval i
        self.ff0_throat_in = np.zeros((int(self.d),int(self.n)))
        self.ff1_throat_in = np.zeros((int(self.d),int(self.n)))
        self.ff2_throat_in = np.zeros((int(self.d),int(self.n)))


    def Calc_Noz_Area(self):
        #CALC_NOZ_AREA determines nozzle area at grid-points x
        self.A = np.zeros((len(self.x)))

        # Nozzle throat located at x = 1.5 in formula below
        for i in range(len(self.A)):
            self.A[i] = 1 + 2.2*( self.x[i] - 1.5 )**(2)


    def Calc_ExactResultsmSW(self):
        # assign values for useful parameters

        fac1 = (self.Gamma - 1)/2

        ithroat = (self.Tot_X_Pts + 1)/2   # grid-point index for nozzle throat
        istart = ithroat + 1   # grid-point index for start of supersonic region
        istop = self.Tot_X_Pts       # grid-point index for end of supersonic region

        # will use expon as a loop parameter so must be an integer

        expon = int((self.Gamma + 1)/(self.Gamma - 1))

        expo2 = -1/(self.Gamma - 1)
        expo1 = self.Gamma*expo2

        # Initialize array to store return arrays
        self.Mach_E = np.zeros((self.Tot_X_Pts))
        self.Mrho_E = np.zeros((self.Tot_X_Pts))
        self.Press_E = np.zeros((self.Tot_X_Pts))
        self.Temp_E = np.zeros((self.Tot_X_Pts))
        self.Vel_E = np.zeros((self.Tot_X_Pts))

        # throat is always at Mach 1 so ...

        self.Mach_E[ithroat] = 1

        # evaluate other flow variables at throat

        argum = 1 + fac1*(self.Mach_E[ithroat])**(2)

        self.Press_E[ithroat] = (argum)**(expo1)
        self.Mrho_E[ithroat] = (argum)**(expo2)
        self.Temp_E[ithroat] = 1/argum

        self.Vel_E[ithroat] = self.Mach_E[ithroat]*np.sqrt(self.Temp_E[ithroat])

        # begin loop over grid-points. Note since nozzle area is symmetric about
        #   the nozzle throat only need to loop over the supersonic flow region

        for i in range(istart, istop + 1):
            # assign value to coefficient of z in polynomial resulting from
            #  area-Mach number (A-M) relation

            c1 = 18750 - 46656*(self.A[i])**(2)

            # specify A-M relation polynomial

            p = [1, 30, 375, 2500, 9375, c1, 15625]

            # evaluate roots of p

            z = np.roots(p)

            # loop over roots: if a root is real and greater than 1 then assign
            #      its square root to M(i) if root real and less than 1, assign
            #      its square root to M(isubson), where isubson = 2*ithroat - i

            for j in range(expon):
                if np.isreal(z[j]):
                    if z(j) > 1:
                        self.Mach_E[i] = np.sqrt(z[j])

                        # calculate other flow variables at i

                        argum = 1 + fac1*z[j]

                        self.Press_E[i] = (argum)**(expo1)
                        self.Mrho_E[i] = (argum)**(expo2)
                        self.Temp_E[i] = 1/argum
                        self.Vel_E[i] = self.Mach_E[i]*np.sqrt(self.Temp_E[i])
                    elif z[j] < 1:
                        isubson = 2*ithroat - i

                        self.Mach_E[isubson] = np.sqrt(z[j])

                        argum = 1 + fac1*z[j]

                        self.Press_E[isubson] = (argum)**(expo1)
                        self.Mrho_E[isubson] = (argum)**(expo2)
                        self.Temp_E[isubson] = 1/argum
                        self.Vel_E[isubson] = self.Mach_E[isubson]*np.sqrt(self.Temp_E[isubson])


    def Calc_ExactResultspSW(self):
        #CALC_EXACTRESULTSPSW assigns exact flow values when shock-wave present

        # Initialize array to store return arrays
        self.Mach_E = np.zeros((self.Tot_X_Pts))
        self.Mrho_E = np.zeros((self.Tot_X_Pts))
        self.Press_E = np.zeros((self.Tot_X_Pts))
        self.Temp_E = np.zeros((self.Tot_X_Pts))
        self.Vel_E = np.zeros((self.Tot_X_Pts))

        # assign values for useful parameters
        ithroat = int((self.Tot_X_Pts + 1)/2)   # grid-point index for nozzle throat
        itp1 = ithroat + 1   # grid-point index for start of supersonic region
        itm1 = ithroat - 1       # grid-point index for end of subsonic region

        fac1 = -1/(self.Gamma -1)
        fac2 = (fac1)**(2)
        fac3 = 2/(self.Gamma - 1)
        fac4 = 2/(self.Gamma + 1)
        fac5 = (self.Gamma - 1)/2

        exponent1 = (self.Gamma + 1)/(self.Gamma - 1)
        exponent2 = self.Gamma/(self.Gamma - 1)
        exponent3 = exponent1/2
        exponent4 = 2*( 2*exponent2 - 1 )

        ratioPA = self.Exit_Pressure*self.A[self.Tot_X_Pts - 1] # Eq.(7.126a) in CFD

        arg1 = fac2 + fac3*( fac4**(exponent1) )*(1/ratioPA)**(2)

        self.Mach_E[self.Tot_X_Pts - 1] = np.sqrt( fac1 + np.sqrt(arg1) ) # Eq.(5.28) in MCF

        arg2 = 1 + fac5*( self.Mach_E[self.Tot_X_Pts - 1] )**(2)

        ratioPoe = (arg2)**(exponent2)      # Eq.(7.128) in CFD
        ratioP21 = ratioPoe*self.Exit_Pressure  # Eq.(7.132) in CFD

        # throat is always at Mach 1 so ...
        self.Mach_E[ithroat - 1] = 1        # Eq.(5.20) in MCF

        # evaluate other flow variables at throat
        argum = 1 + fac5*(self.Mach_E[ithroat- 1])**(2)

        self.Press_E[ithroat - 1] = (argum)**(-1 * exponent2)     # Eq.(3.30) in MCF
        self.Mrho_E[ithroat - 1] = (argum)**(fac1)            # Eq.(3.31) in MCF
        self.Temp_E[ithroat - 1] = 1/argum                   # Eq.(3.28) in MCF

        self.Vel_E[ithroat - 1] = self.Mach_E[ithroat - 1]*np.sqrt(self.Temp_E[ithroat - 1]) # Eq. (7.76) in CFD

        # begin calculation to locate shock-wave position
        term1 = (6**(exponent4))*(1/ratioP21)**(2)

        fcoeff7 = 855100224 - term1

        # assign polynomial that yields Mach number just before the shock-wave 
        #   see Eq.(7.127) in CFD. NOTE: flow is supersonic there!
        f = [ 16807, 576240, 8406930, 67347560, 316914885, fcoeff7, 1123343340, 214527600, -1 * 594129375, 216650000, -1 * 34518750, 2625000, -1 * 78125 ]

        z = np.roots(f)

        # loop over roots: if root is real, and greater than 1 (supersonic flwo), 
        #   then assign it to z1 and set M1 = sqrt(z1)
        for j in range(int(exponent4)):
            if np.isreal(z[j]):
                if z[j] > 1:
                    z1 = z[j]
                    M1 = np.sqrt(z1)

        # calculate jump in pressure and Mach number across shock-wave
        self.SW_JumpP = 1 + ((2*self.Gamma)/(self.Gamma+1))*( M1**2 - 1 )  # Eq.(7.137) in CFD
        self.SW_JumpM = np.sqrt( (1 + fac5*M1**2)/(self.Gamma*M1**2 - fac5) ) # Eq.(7.138) CFD

        arg3 = 1 + fac5*z1
        arg4 = fac4*arg3

        # calculate nozzle area Asw at shock-wave location
        Asw = (1/M1)*( arg4 )**(exponent3)     # Eq.(5.20) in MCF

        # calculate position Xsw of shock-wave - use root greater than 1.5 since
        #   shock-wave is in divergent part of nozzle
        Xsw = 1.5 + np.sqrt( (Asw - 1)/2.2 )     # solve Eq.(7.73) in CFD for x

        # calculate grid-point index isw of shock-wave location
        isw = np.round( 1 + Xsw/self.Del_x)

        # define grid-point label for point just downstream of shock-wave
        iswp1 = isw + 1

        # now that we know where the shock-wave is we can determine the exact
        #   values of the flow variables at all grid-points.
        # begin by calculating flow variables upstream of shock-wave. First
        #   need to determine Mach number at all upstream grid-points
        # flow is subsonic before throat
        for i in range(itm1):
            term2 = (6**(exponent1))*(self.A[i]**(2))

            gcoeff1 = 18750 - term2

            g = [ 1, 30, 375, 2500, 9375, gcoeff1, 15625 ]

            z = np.roots(g)

            # want root that's real, positive, and less than 1 (subsonic flow)
            for j in range(int(exponent1)):
                if np.isreal(z[j]):
                    if z[j] > 0:
                        if z[j] < 1:
                            self.Mach_E[i] = np.sqrt(z[j])

        # now repeat for supersonic flow upstream of shock-wave
        for i in range(int(itp1) - 1, int(isw)):
            term2 = ( 6**(exponent1) )*( self.A[i]**(2) )

            gcoeff1 = 18750 - term2

            g = [ 1, 30, 375, 2500, 9375, gcoeff1, 15625 ]

            z = np.roots(g)

            # want root that's real and greater than 1 (supersonic flow)
            for j in range(int(exponent1)):
                if np.isreal(z[j]):
                    if z[j] > 1:
                        self.Mach_E[i] = np.sqrt(z[j])

        # now calculate flow variables upstream of shock-wave
        for i in range(int(isw)):
            arg5 = 1 + fac5*(self.Mach_E[i])**(2)

            self.Press_E[i] = (arg5)**(-1 * exponent2)
            self.Mrho_E[i] = (arg5)**(fac1)
            self.Temp_E[i] = 1/arg5

            self.Vel_E[i] = self.Mach_E[i]*np.sqrt(self.Temp_E[i])


        # now calculate flow downstream of shock-wave flow is subsonic here
        # start by calculating sonic area AstarE for nozzle exit
        arg7 = fac4*( 1 + fac5*( self.Mach_E[self.Tot_X_Pts - 1] )**(2) )

        # derivation of following formula in my research notebook, 07-31-19 entry,
        #   page 11
        AstarE = self.A[self.Tot_X_Pts - 1]*self.Mach_E[self.Tot_X_Pts - 1]*( arg7 )**(-1 * exponent3)

        # begin loop over downstream grid-points
        for i in range(int(iswp1)-1, self.Tot_X_Pts):
            term3 = ( (6**(exponent1))*(self.A[i])**(2) )/( AstarE**(2) )

            hcoeff1 = 18750 - term3

            h = [1, 30, 375, 2500, 9375, hcoeff1, 15625 ]

            z = np.roots(h)

            # want root that's real, positive , and less than 1 (subsonic flow)
            for j in range(int(exponent1)):
                if np.isreal(z[j]):
                    if z[j] > 0:
                        if z[j] < 1:
                            self.Mach_E[i] = np.sqrt(z[j])


        # calculate flow variables downstream of shock-wave
        for i in range(int(iswp1) -1, self.Tot_X_Pts):
            arg6 = 1 + fac5*(self.Mach_E[i])**(2)

            # following formulas derived in my research notebook, 07-30-19 entry,
            #      page 7, Eqs.(15a-c)
            self.Press_E[i] = ratioP21*(arg6)**(-1 * exponent2)
            self.Mrho_E[i] = ratioP21*(arg6)**(fac1)
            self.Temp_E[i] = 1/arg6

            self.Vel_E[i] = self.Mach_E[i]*np.sqrt(self.Temp_E[i])


    def SetInCond(self):
        # determine coefficients and set Tiny to small value for later use
        Coef1 = 1 /(self.Gamma -1)
        Coef2 = self.Gamma/2

        ithroat = (self.Tot_X_Pts + 1)/2

        # initialize arrays
        Mrho = np.zeros((self.Tot_X_Pts))    # will store shifted mass density
        Temp = np.zeros((self.Tot_X_Pts))    # will store shifted temperature
        self.InitVal = np.zeros((self.d,self.Tot_X_Pts))       # will store flow variables
        V = np.zeros((self.Tot_X_Pts))       # will store velocity
        Loc_TSteps = np.zeros((self.Tot_X_Pts- 2))    # will store local time-step

        # add random shift to In_Mass_Flow_E
        self.In_Mass_Flow_Noisy = self.In_Mass_Flow*(1 + (-1 * self.ICMFlowErrScale +2*self.ICMFlowErrScale*random.random()))

        # begin loop over grid-points - initialize mass density and temperature
        # Introduce a small shift
        #Tiny = 10**(-8)
        #print("Shock flag =", Shock_Flag)
        for i in range(int(self.Tot_X_Pts)):
            if self.Shock_Flag == 0:
                if ((i == 0) or (i == ithroat - 1)):
                    Mrho[i] = self.Mrho_E[i]
                    Temp[i] = self.Temp_E[i]
        
                elif ((i != 0) and (i != ithroat - 1)):
                    Mrho[i] = self.Mrho_E[i]*( 1 + (-1 *self.ICrhoErrScale + 2*self.ICrhoErrScale*random.random()) )

                    if Mrho[i] > 1:
                        Mrho[i] = 1

                    Temp[i] = self.Temp_E[i]*( 1 + (-1 * self.ICtempErrScale + 2*self.ICtempErrScale*random.random()) )

                    if Temp[i] > 1:
                        Temp[i] = 1

            elif self.Shock_Flag == 1:
                if ((i == 0) or (i == ithroat - 1)):
                    Mrho[i] = self.Mrho_E[i]
                    Temp[i] = self.Temp_E[i]
        
                elif ( (i != 0) and (i != ithroat - 1) ):
                    Mrho[i] = self.Mrho_E[i]*( 1 + (-1 * self.ICrhoErrScale + 2*self.ICrhoErrScale*random.random()) )

                    if Mrho[i] > 1:
                        Mrho[i] = 1

                    Temp[i] = self.Temp_E[i]*( 1 + (-1 * self.ICtempErrScale + 2*self.ICtempErrScale*random.random()) )

                    if Temp[i] > 1:
                        Temp[i] = 1
                
            # assign initial condition to flow variables
            self.InitVal[0, i] = Mrho[i]*self.A[i]
            self.InitVal[1, i] = self.In_Mass_Flow
            self.InitVal[2, i] = self.InitVal[0, i]*( Coef1*Temp[i] + Coef2*(self.InitVal[1, i]/self.InitVal[0, i])**(2) )

        # determine time-step Del_t using Courant-Friedrichs-Levy (CFL) stability
        #   condition with C = 0.5
        C = 0.5

        for i in range(1, (self.Tot_X_Pts - 1)):
            V[i] = self.InitVal[1, i]/self.InitVal[0, i]
            Loc_TSteps[i - 2] = (C*self.Del_x)/(np.sqrt(Temp[i]) + V[i])

        # Del_t is the min value in Loc_TSteps
        self.Delta_t = np.min(Loc_TSteps)


    def InitParms(self):
        #INITPARMS initializes parameters for quantum ODE algorithm
        self.k = 1+ np.ceil(np.log(1/self.err1)/np.log(self.in_n))

        self.N = self.in_n**(self.k-1)

        Tylr_Steps = self.in_n * self.N

        ratio = self.Tot_TSteps/Tylr_Steps

        while ratio > 1:
            self.in_n = self.in_n + 1

            self.k = 1 + np.ceil(np.log(1/self.err1)/np.log(self.in_n))
            self.N = self.in_n**(self.k-1)

            Tylr_Steps = self.in_n * self.N
            ratio = self.Tot_TSteps/Tylr_Steps

        temp_n = (self.Tot_TSteps/ratio)**(1/self.k)

        self.final_n = np.ceil(temp_n)

        self.k = 1 + np.ceil(np.log(1/self.err1)/np.log(self.final_n))

        self.N = (self.final_n)**(self.k-1)

        self.delta1 = 1 - (1-self.delta)**(1/self.N)


    def IPrtn(self, b):# Init, Fnl, nn, NN ): a, b, n, N
        #IPRTN partitions interval [Init, Fnl] into nn*NN subintervals
        hh = (b - self.a)/self.n     # width of a primary subinterval
        self.hbar = hh/self.N             # width of a sub-subinterval

        self.t = np.zeros((int(self.N),int(self.n)))        # create/initialize times array first (second)
        
        # slot specifies sub- (primary) subinterval
        # fill in t array-element values Note t(N,n) = Fnl
        for i in range(int(self.n)):
            for j in range(int(self.N)):
                if i == 0:
                    self.t[j, i] = self.a + j*self.hbar    # fill times in first subinterval
                else:
                    self.t[j, i] = self.t[int(self.N) - 1, i-1] + j*self.hbar # fill remaining subintervals


    def IntegrateODE(self):#d, n, N, hbar, r, Del_x, Gamma, Tot_Int_Pts, k, Tot_X_Pts, Shock_Flag, Exit_Pressure, ithroat, a, delta1, rho, InitVal, A, t, U2_in, ff0_throat_in, ff1_throat_in, ff2_throat_in, Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, In_Mass_Flow):
        #INTEGRATEODE numerically integrates ODE for 1D Navier-Stokes flow

        # initialize arrays ff0_throat, ff1_throat, ff2_throat to values passed
        #   through function input arguments 
        ff0_throat = self.ff0_throat_in
        ff1_throat = self.ff1_throat_in
        ff2_throat = self.ff2_throat_in

        # similarly, initialize U2 to U2_in
        U2 = self.U2_in

        # deleting things that are no longer needed
        del self.ff0_throat_in, self.ff1_throat_in, self.ff2_throat_in, self.U2_in

        # Begin loop over the subintervals i result is approximate solution z(t).
        for i in range(int(self.n)):
            #build Taylor polynomials l**{s}_{i}(t)for subinterval i at all
            #   interior grid-points store polynomial coefficients in 
            #       StoreLz(d,r+2,Tot_Int_Pts,N)
            #   ff_throat is a d x (r+1) array storing the values of residual and
            #       its first r time derivatives at the nozzle throat at the end 
            #       of subinterval i
            StoreLz, ff_throat = BldTPoly(d , n, N, hbar, r, InitVal, Del_x, Gamma, Tot_Int_Pts,Tot_X_Pts, A, Shock_Flag,Exit_Pressure, ithroat)

            # store values of ff_throat for subinterval i for easy of writing to
            #   files

            ff0_throat[:, i] = np.abs(ff_throat[:, 0])
            ff1_throat[:, i] = np.abs(ff_throat[:, 1])
            ff2_throat[:, i] = np.abs(ff_throat[:, 2])

            # store N intermediate times for subinterval i

            StoreTimes4i = t[:, i]

            # define Start time for subinterval i

            if i == 0:
                Start = a
            else:
                Start = t[N - 1, i-1]


            #  gInt = d x Tot_Int_Pts array storing integral of each component of
            #               g_ij over subinterval i for each interior grid-point
            #

            gInt = IntegrateGij(StoreLz, StoreTimes4i, Start, d, r, N, delta1, hbar, rho, Tot_Int_Pts, A, Gamma, Del_x, Shock_Flag, i, Exit_Pressure)

            # add gInt for subint i to InitVal to get InitVal for 
            #   subinterval i + 1 at each interior grid-point

            for ll in range(1, (Tot_Int_Pts+1)):
                IP_Label = ll - 1

                for m in range(d):
                    InitVal[m, ll] = InitVal[m, ll] + gInt[m, IP_Label]


            # use flow boundary conditions to determine InitVal for
            #   subinterval i + 1 at two boundary points 1 and Tot_X_Pts. Store
            #   them in InitVal_BVals

            if Shock_Flag == 0:
                InitVal_BVals = CalcBCmSW(InitVal,A,Gamma,d,Tot_X_Pts)
            elif Shock_Flag == 1:
                InitVal_BVals = CalcBCpSW(InitVal,A,Gamma,d,Tot_X_Pts, Exit_Pressure)
            else:
                print('Unknown Shock_Flag value: ', Shock_Flag)


            # transfer boundary values from InitVal_Bvals to InitVal. Resulting 
            #       InitVal then contains initial values of primary flow 
            #       variable U for subinterval i + 1 at each grid-point.

            for m in range(d):
                InitVal[m, 0] = InitVal_BVals[m, 0]
                InitVal[m, Tot_X_Pts] = InitVal_BVals[m, 1]


            # store U2 values at start of subinterval i + 1 at all gridpoints

            for gridpt in range(Tot_X_Pts):
                U2[gridpt, i+1] = InitVal[1, gridpt]


            # evaluate physical flow variables from InitVal for next subinterval

            Mach_D, Mrho_D, Press_D, Temp_D, Vel_D = Calc_FlowVarResults(Gamma,Tot_X_Pts,A,InitVal)

            # output initial simulation flow variables for subinterval i+1

            print('Code has completed subinterval ',i)

            if i != n:
                print('  Initial condition for next subinterval is:')
            elif i == n:
                print('  Final result for steady state U values are:')


            print(np.transpose(InitVal))

            nextStartTime = n**(k-1)*hbar*i

            print('Next subint start-time = ', nextStartTime)


        # calculate relative error in primary flow variables
        self.Rel_MachErr = np.divide(np.abs(Mach_D - Mach_E), Mach_E)
        self.Rel_MrhoErr = np.divide(np.abs(Mrho_D - Mrho_E), Mrho_E)
        self.Rel_PressErr = np.divide(np.abs(Press_D - Press_E), Press_E)
        self.Rel_TempErr = np.divide(np.abs(Temp_D - Temp_E), Temp_E)
        self.Rel_VelErr = np.divide(np.abs(Vel_D - Vel_E), Vel_E)

        # calculate mean mass flow rate at final time define array to store
        self.MeanU2 = np.mean(U2[:, n+1])
        self.AvU2 = np.zeros(Tot_X_Pts)

        # calculate mean and standard deviation of relative errors and store
        self.MeanRelTempErr = np.mean(self.Rel_TempErr)
        self.MeanRelMachErr = np.mean(self.Rel_MachErr)
        self.MeanRelMrhoErr = np.mean(self.Rel_MrhoErr)
        self.MeanRelPressErr = np.mean(self.Rel_PressErr)

        self.SDevRelTempErr = np.std(self.Rel_TempErr)
        self.SDevRelMachErr = np.std(self.Rel_MachErr)
        self.SDevRelMrhoErr = np.std(self.Rel_MrhoErr)
        self.SDevRelPressErr = np.std(self.Rel_PressErr)

        self.AvRelTempErr = np.zeros(Tot_X_Pts)
        self.AvRelMachErr = np.zeros(Tot_X_Pts)
        self.AvRelMrhoErr = np.zeros(Tot_X_Pts)
        self.AvRelPressErr = np.zeros(Tot_X_Pts)

        self.AvPlusSDevRelTempErr = np.zeros(Tot_X_Pts)
        self.AvPlusSDevRelMachErr = np.zeros(Tot_X_Pts)
        self.AvPlusSDevRelMrhoErr = np.zeros(Tot_X_Pts)
        self.AvPlusSDevRelPressErr = np.zeros(Tot_X_Pts)

        self.AvMinusSDevRelTempErr = np.zeros(Tot_X_Pts)
        self.AvMinusSDevRelMachErr = np.zeros(Tot_X_Pts)
        self.AvMinusSDevRelMrhoErr = np.zeros(Tot_X_Pts)
        self.AvMinusSDevRelPressErr = np.zeros(Tot_X_Pts)

        for col in range(Tot_X_Pts):
            self.AvU2[col] = self.MeanU2

            self.AvRelTempErr[col] = self.MeanRelTempErr
            self.AvRelMachErr[col] = self.MeanRelMachErr
            self.AvRelMrhoErr[col] = self.MeanRelMrhoErr
            self.AvRelPressErr[col] = self.MeanRelPressErr

            self.AvPlusSDevRelTempErr[col] = self.MeanRelTempErr + self.SDevRelTempErr
            self.AvPlusSDevRelMachErr[col] = self.MeanRelMachErr + self.SDevRelMachErr
            self.AvPlusSDevRelMrhoErr[col] = self.MeanRelMrhoErr + self.SDevRelMrhoErr
            self.AvPlusSDevRelPressErr[col] = self.MeanRelPressErr + self.SDevRelPressErr
                                
            self.AvMinusSDevRelTempErr[col] = self.MeanRelTempErr - self.SDevRelTempErr
            self.AvMinusSDevRelMachErr[col] = self.MeanRelMachErr - self.SDevRelMachErr
            self.AvMinusSDevRelMrhoErr[col] = self.MeanRelMrhoErr - self.SDevRelMrhoErr
            self.AvMinusSDevRelPressErr[col] = self.MeanRelPressErr - self.SDevRelPressErr


        #return U2, Mach_D, Mrho_D, Press_D, Temp_D, Vel_D, ff0_throat, ff1_throat, ff2_throat


    def BldTPoly(self):#dd,nn,NN,hb,rr,InVal,Del_x, Gamma,Tot_Int_Pts, Tot_X_Pts, A, Shock_Flag, Exit_Pressure, ithroat):
        #BLDTPOLY calculates all Taylor polynomial coefficients for subint i
        
        # ll will store rmaxp1 Taylor poly. coeffs. for each component of ODE 
        # driver function, subsubinterval j, and interior grid-point.
        ll = np.zeros((int(self.d), int(self.r + 2), int(self.Tot_Int_Pts), int(self.N)))   
        
        for j in range(int(self.N)):
            # evaluate ODE driver func. and first rr time derivatives at InVal
            #    store in ff = dd x (rr+1) x Tot_Int_Pts array
            ff = self.Derivs()

            # for each subsubinterval j, calculate and store Taylor poly coeffs. 
            #  for each component of ODE driver function and interior grid-point  
            #  in array mat = dd x (rr+2) x Tot_Int_Pts 
            mat = BldTMat(ff,dd,rr,InVal,Tot_Int_Pts,Tot_X_Pts)  

            ll[:, :, :, int(j)] = mat   # transfer Taylor polys coeffs for each component
                        # of driver function, subsubint j and interior
                        # grid-point

            
            InVal = NextInCond(mat,InVal,hb,dd,rr,A,Gamma,Tot_Int_Pts, Tot_X_Pts,Shock_Flag,Exit_Pressure) 

        # if change store values of ff at ithroat to send back to 
        #   main program on last iteration
        return ll, ff[:, :, int(self.ithroat) - 1]

    
    def Derivs(self):# d,r,InitVal,Del_x,Gamma,Tot_Int_Pts,Tot_X_Pts, A,Shock_Flag):
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
        F = self.CalcFlux()


        # 1.b Evaluate source current J in momentum equation of motion which is
        #       a 1 x Tot_Int_Pts array
        J = self.CalcSource()

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
        if self.Shock_Flag == 0: ff_Bvals = self.CalcfBvalsmSW(ff_vals)
        elif self.Shock_Flag == 1: ff_Bvals = self.CalcfBvalspSW(ff_vals)
        else: print('Unkown Shock_Flag value: ', self.Shock_Flag)

        #   2.b calculate first time derivative of flow fluxes at all grid-points
        #       
        #       dFdt = d x Tot_X_Pts array storing the flow flux first time
        #                               derivatives at all grid-points.
        dFdt = self.Calc_dFdt(ff_vals, ff_Bvals)

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
        if self.Shock_Flag == 0: dffdt_Bvals = CalcdfdtBvalsmSW(InitVal,Gamma,ff_Bvals,dffdt_vals,d,Tot_Int_Pts)
        elif self.Shock_Flag == 1: dffdt_Bvals = CalcdfdtBvalspSW(InitVal,Gamma,ff_Bvals,dffdt_vals,d, Tot_Int_Pts)
        else: print('Unknown Shock_Flag value: ', self.Shock_Flag)
        

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
        for ll in range(self.Tot_Int_Pts):
            for k in range(self.d):
                ff[k, 2, ll] = d2ffdt2_vals[k, ll]

        return ff


    def CalcFlux(self):
        #CALCFLUX evaluates flow fluxes at all grid-points for Nav-Stokes dynamics
        F = np.zeros((self.d, self.Tot_X_Pts))

        fac1 = (3 - self.Gamma)/2
        fac2 = (self.Gamma - 1)/self.Gamma
        fac3 = (self.Gamma - 1)/2

        # calculate fluxes loop over all grid-points 
        for i in range(self.Tot_X_Pts):
            F[0, i] = self.InitVal[1, i]

            # introduce some useful definitions
            ratio1 = 1/self.InitVal[0, i]
            ratio2 = (self.Gamma*self.InitVal[1, i])/(self.InitVal[0, i]*self.InitVal[0, i])

            term1 = self.InitVal[0, i]*self.InitVal[2, i]
            term2 = self.InitVal[1, i]*self.InitVal[1, i]

            F[1, i] = ratio1*( fac2*term1 + fac1*term2 )
            F[2, i] = ratio2*( term1 - fac3*term2 )

        return F


    def CalcSource(self):#U, A, Gamma, Tot_X_Pts):
        #CALCSOURCE evaluates the source current in the Navier-Stokes dynamics
        J = np.zeros(self.Tot_Int_Pts)

        # calculate source current loop over interior grid-points with
        #       2 <= ll <= Tot_Int_Pts +1
        for i in range(1, (self.Tot_Int_Pts + 1)):    # ll indexes the grid-points                                       
            fac1 = ((self.Gamma - 1)/self.Gamma) * np.log(self.A[i+1]/self.A[i-1])
            fac2 = self.Gamma/2

            term1 = self.InitVal[0, i] * self.InitVal[2, i]
            term2 = self.InitVal[1, i] * self.InitVal[1, i]
                                    
            ratio1 = 1/self.InitVal[0, i]

            J[i - 1] = fac1*ratio1*( term1 - fac2*term2 )


        return J


    def CalcFunc(self, F, J):#, Del_x, d, Tot_Int_Pts ):
        #CALCFUNC evaluates the ODE driver function for 1D Nav.-Stokes dynamics
        ff_vals = np.zeros((self.d, self.Tot_Int_Pts))

        Fac1 = (-1)/(2*self.Del_x)

        # evaluate ODE driver function, loop over interior grid-points
        #   NOTE: 1. gridLabel gives grid-point label for each interior grid-point
        #         2. Since F is specified at all grid-points we must use gridLabel
        #               to get the values of F to left and right of an interior
        #               grid-point
        #         3. J is specified only at interior points so i is used to get
        #               values of j at specified interior point
        for i in range(self.Tot_Int_Pts):
            gridLabel = i + 1

            ff_vals[0, i] = Fac1*( F[0, gridLabel + 1] - F[0, gridLabel - 1])
            ff_vals[1, i] = Fac1*( (F[1, gridLabel+1] - F[1, gridLabel-1]) - J[i] )
            ff_vals[2, i] = Fac1*( F[2, gridLabel + 1] - F[2, gridLabel - 1] )


        return ff_vals


    def CalcfBvalsmSW(self, ff_vals):#U,Gamma,ff_vals,d,Tot_Int_Pts):
        #CALCFBVALSMSW assigns ODE driver function boundary values - no shock-wave
        ff_Bvals = np.zeros((self.d,2))

        Fac3 = self.InitVal[1, 0]/self.InitVal[0, 0]

        # calculate ff_Bvals at nozzle entrance
        ff_Bvals[0, 0] = 0
        ff_Bvals[1, 0] = 2*ff_vals[1, 0] - ff_vals[1, 1]
        ff_Bvals[2, 0] = self.Gamma*Fac3*ff_Bvals[1, 0]

        # calculate use boundary conditions on flow variables to determine
        #   ff_Bvals at nozzle exit
        for k in range(3): ff_Bvals[k, 1] = 2*ff_vals[k, self.Tot_Int_Pts] - ff_vals[k, (self.Tot_Int_Pts-1)]

        return ff_Bvals


    def CalcfBvalspSW(self, ff_vals):
        #CALCFBVALSPSW assigns boundary values to ODE driver fcn - with shock-wave 
        ff_Bvals = np.zeros((self.d,2))

        Fac1 = self.Gamma/2
        Fac3 = self.InitVal[1, 0]/self.InitVal[0, 0]

        # calculate ff_Bvals at nozzle entrance
        ff_Bvals[0, 0] = 0
        ff_Bvals[1, 0] = 2*ff_vals[1, 0] - ff_vals[1, 1]
        ff_Bvals[2, 0] = self.Gamma*Fac3*ff_Bvals[1, 0]

        # calculate use boundary conditions on flow variables to determine
        #   ff_Bvals at nozzle exit
        for k in range(self.d):
            ff_Bvals[k, 1] = 2*ff_vals[k, self.Tot_Int_Pts - 1] - ff_vals[k, self.Tot_Int_Pts-2]

        ff_Bvals[2, 1] = Fac1*Fac3*( 2*ff_Bvals[1, 1] - Fac3*ff_Bvals[0, 1] )

        return ff_Bvals


    def Calc_dFdt(self, ff_vals, ff_Bvals):
        #CALC_DFDT evaluates first time derivative of flow fluxes at grid-points
        dFdt = np.zeros((self.d, self.Tot_X_Pts))

        TotXPtm1 = self.Tot_X_Pts - 1

        fac4 = self.Gamma -1
        fac5 = self.Gamma/2
        fac6 = fac4*fac5
        fac7 = (3-self.Gamma)/2
        fac8 = fac4/self.Gamma

        # evaluate dFdt at boundaries:
        for ll in range(2):
            if ll == 0: # at nozzle entrance
                fac1 = self.InitVal[1, 0]/self.InitVal[0, 0]
                fac2 = self.InitVal[2, 0]/self.InitVal[0, 0]
                fac3 = fac1*fac2
                fac9 = (fac1)**(2)

                dFdt[0, ll] = ff_Bvals[2, ll]

                dFdt[1, ll] = fac7*fac1*(2*ff_Bvals[1, ll]-fac1*ff_Bvals[0, ll]) + fac8*ff_Bvals[2, ll]

                dFdt[2, ll] = self.Gamma*(fac1*ff_Bvals[2, ll]+fac2*ff_Bvals[1, ll]-fac3*ff_Bvals[0, ll])-fac6*fac9*(3*ff_Bvals[1, ll])-2*fac1*ff_Bvals[0, ll]
            elif ll == 1: # at nozzle exit
                fac1 = self.InitVal[1, self.Tot_X_Pts - 1]/self.InitVal[0, self.Tot_X_Pts - 1]
                fac2 = self.InitVal[2, self.Tot_X_Pts - 1]/self.InitVal[0, self.Tot_X_Pts - 1]
                fac3 = fac1*fac2
                fac9 = (fac1)**(2)

                dFdt[0, self.Tot_X_Pts- 1] = ff_Bvals[1, ll]

                dFdt[1, self.Tot_X_Pts- 1]= fac7*fac1*(2*ff_Bvals[1, ll]-fac1*ff_Bvals[0, ll]) + fac8*ff_Bvals[2, ll]

                dFdt[2, self.Tot_X_Pts- 1] = self.Gamma*(fac1*ff_Bvals[2, ll] + fac2*ff_Bvals[2, ll]- fac3*ff_Bvals[0, ll])-fac6*fac9*(3*ff_Bvals[1, ll]-2*fac1*ff_Bvals[0, ll])
            else:
                print('Unknown switch case label: ', ll)

        # evaluate dFdt at interior grid-points
        for ll in range(1, TotXPtm1 + 1):
            fac1 = self.InitVal[1, ll]/self.InitVal[0, ll]
            fac2 = self.InitVal[2, ll]/self.InitVal[0, ll]
            fac3 = fac1*fac2
            fac9 = (fac1)**(2)

            IPLabel = ll -1

            dFdt[0, ll] = ff_vals[1, IPLabel- 1]

            dFdt[1, ll] = fac7*fac1*(2*ff_vals[1, IPLabel- 1]-fac1*ff_vals[0, IPLabel- 1])+ fac8*ff_vals[2, IPLabel- 1]

            dFdt[2, ll] = self.Gamma*(fac1*ff_vals[2, IPLabel- 1] + fac2*ff_vals[1, IPLabel- 1] - fac3*ff_vals[1, IPLabel- 1])-fac6*fac9*(3*ff_vals[1, IPLabel- 1]-2*fac1*ff_vals[0, IPLabel- 1])

        return dFdt


inst = ns_q()