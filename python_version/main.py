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
        for item in dir(self):
            if "__" not in item and not callable(eval("self." +item)):
                if val:
                    print(item, "=", eval("self." + item))
                else:
                    print(item)


    def __init__(self):
        self.InitCalcParms()

        start = time.time()

        # integrate ODE
        '''
        U2, Mach_D, Mrho_D, Press_D, Temp_D, Vel_D, Rel_MachErr, Rel_MrhoErr, Rel_PressErr, Rel_TempErr, Rel_VelErr, AvRelTempErr, AvPlusSDevRelTempErr, AvMinusSDevRelTempErr, AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr, AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr, AvRelPressErr, AvPlusSDevRelPressErr, AvMinusSDevRelPressErr, AvU2, ff0_throat, ff1_throat, ff2_throat = IntegrateODE(d, n, N, hbar, r, Del_x, Gamma, Tot_Int_Pts, k, Tot_X_Pts, Shock_Flag, Exit_Pressure, ithroat, a, delta1, rho, InitVal, A, t, U2_in, ff0_throat_in, ff1_throat_in, ff2_throat_in, Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, In_Mass_Flow)

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
        self.check("d, n, N, hbar, r, Del_x, Gamma, Tot_Int_Pts, k, Tot_X_Pts, Shock_Flag, Exit_Pressure, ithroat, a, delta1, rho, InitVal, A, t, U2_in, ff0_throat_in, ff1_throat_in, ff2_throat_in, Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, In_Mass_Flow")

    def InitCalcParms(self):
        #d, n, N, hbar, r, Del_x, Gamma, Tot_Int_Pts, k, Tot_X_Pts, Shock_Flag, Exit_Pressure, ithroat, a, delta1, rho, InitVal, A, t, U2_in, ff0_throat_in, ff1_throat_in, ff2_throat_in, Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, In_Mass_Flow

        # set up x-grid for problem

        self.x = np.linspace(self.x_min, self.x_max, self.Tot_X_Pts)
        self.Del_x = self.x[1] - self.x[0]

        # calculate nozzle area at each grid-point
        self.Calc_Noz_Area()

        # calculate exact steady-state solution of primary flow variables: Mach
        #   number, mass density, pressure, and temperature
        if self.Shock_Flag == 0:
            self.Calc_ExactResultsmSW()
        elif self.Shock_Flag == 1:
            self.Calc_ExactResultspSW()
        else:
            print('Unknown Shock_Flag value: ', self.Shock_Flag)

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

inst = ns_q()