from qiskit.algorithms import AmplitudeEstimation, EstimationProblem
from qiskit.circuit import QuantumCircuit
from qiskit import Aer, IBMQ

import time, warnings, random, math, ctypes, os
from timeit import default_timer as timer
from numba import cuda, jit, njit, prange
from pprint import pprint
import numpy as np

warnings.filterwarnings("ignore")


"""
This to look into:
    - Look at ways to use GijVals in IntergrateGij better, that is where a lot of time is spent

"""


class BernoulliA(QuantumCircuit):
    """
    A circuit representing the Bernoulli A operator.
    
    """
    def __init__(self, probability):
        super().__init__(1)  # circuit on 1 qubit

        theta_p = 2 * np.arcsin(np.sqrt(probability))
        self.ry(theta_p, 0)


class BernoulliQ(QuantumCircuit):
    """
    A circuit representing the Bernoulli Q operator.
    
    """
    def __init__(self, probability):
        super().__init__(1)  # circuit on 1 qubit

        self._theta_p = 2 * np.arcsin(np.sqrt(probability))
        self.ry(2 * self._theta_p, 0)

    def power(self, k):
        # implement the efficient power of Q
        q_k = QuantumCircuit(1)
        q_k.ry(2 * k * self._theta_p, 0)
        return q_k


class ns_q:
    # controls for features of the execution of the code
    log_QAmpEst = True # log data from the quantum algorithm, "True" if you want data to be logged
    run_time_option = 2 # 0 = original method, 1 = qasm simulator, 2 = real machine

    # calculation controls
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


    def __init__(self):
        start = time.time()
        self.InitCalcParms()

        # integrate ODE
        self.IntegrateODE()

        # stop the timer for the calculation
        runtime = (time.time() - start)/60

        print('Program runtime (minutes) = ', runtime)
        print('Program runtime per subinterval(minutes) = ', runtime/self.n)


    def InitCalcParms(self):
        print("Initializing")

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

        # sets the directory where all of the data generated will end up
        self.data_dir = f"{time.ctime().replace(' ', '_')}"

        # if data from the algo should be logged, it is
        if self.log_QAmpEst:
            self.QAmpEst_log_f = os.path.join(self.data_dir, "QAmpEst")
            
            # holds omega and result values
            self.data_vals = []

        if self.run_time_option == 0:
            # loads the dll that runs the randQAEA functions
            self.c_lib = ctypes.CDLL("./c++_backend/libbackend.so").QAmpEst
            self.c_lib.restype = ctypes.c_double

            self.algo_func = self.original_algo

        elif self.run_time_option == 1:
            self.ae = AmplitudeEstimation(
                num_eval_qubits = 4,  # the number of evaluation qubits specifies circuit width and accuracy
                quantum_instance = Aer.get_backend('qasm_simulator')
            )

            self.algo_func = self.simulator_algo
        
        elif self.run_time_option == 2:
            # loading my account
            IBMQ.load_account()

            self.ae = AmplitudeEstimation(
                num_eval_qubits = 4,  # the number of evaluation qubits specifies circuit width and accuracy
                quantum_instance = IBMQ.get_provider('ibm-q').get_backend('ibmq_belem')
            )

            self.algo_func = self.real_algo

        else: print("Invalid Runtime Option")

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
        self.U2 = np.zeros((int(self.Tot_X_Pts), int(self.n+1)))

        # use initial condition to assign U2 at start of subinterval i = 1
        for gridpt in range(int(self.Tot_X_Pts)):
            self.U2[gridpt, 0] = self.InitVal[1, gridpt]

        # define ithroat as grid-point index of nozzle throat
        self.ithroat = (self.Tot_X_Pts + 1)/2

        # calculate total number of amplitude estimates needed (TotRuns)
        FudgeFactor = 1.25
        TempTot = np.ceil(FudgeFactor*(-8 * np.log(self.delta1)))

        if (TempTot % 2 == 0): self.TotRuns = TempTot + 1
        else: self.TotRuns = TempTot

        # used to test for division by zero below
        self.Tolerance = 10**(-1*12)

        self.Mach_D = np.zeros((self.Tot_X_Pts))
        self.Mrho_D = np.zeros((self.Tot_X_Pts))
        self.Press_D = np.zeros((self.Tot_X_Pts))
        self.Temp_D = np.zeros((self.Tot_X_Pts))
        self.Vel_D = np.zeros((self.Tot_X_Pts))

        # initialize arrays to store value of driver function f(z(t)) and its first
        #   r time derivatives at nozzle throat at end of each subinterval i
        self.ff0_throat = np.zeros((int(self.d),int(self.n)))
        self.ff1_throat = np.zeros((int(self.d),int(self.n)))
        self.ff2_throat = np.zeros((int(self.d),int(self.n)))

        # nice messages to the user
        print(np.transpose(self.InitVal))
        print('Above are initial values of computational flow values U')
        input("\n--> Press enter to see more runtime values:")
        print()        

        print("Run time method")
        if self.run_time_option == 0: print(" -- run time option = original")
        elif self.run_time_option == 1: print("-- run time option = qasm simulator")
        else: print("-- run time option = real machine")
        print()

        print("Algorithm parameters")
        print('-- Delta_t =', self.Delta_t)
        print('-- hbar =', self.hbar)
        print('-- b =', b)
        print('-- n =', self.n)
        print('-- k =', self.k)
        print('-- In_Mass_Flow_Noisy =', self.In_Mass_Flow_Noisy)
        
        if self.Shock_Flag == 1:
            print('-- SW_JumpP =', self.SW_JumpP)
            print('-- SW_JumpM =', self.SW_JumpM)
        
        print()
        print("Data handling values")
        print(f"-- Log data from QAmpEst = {self.log_QAmpEst}", end="")
        if self.log_QAmpEst: print(f", logged in --> {self.data_dir}", end="")
        print()

        input("\n--> Press enter to continue run the algorithm:")
        print()

        os.mkdir(self.data_dir) # makes the data directory
        self.N = int(self.N)

        # has to be here
        if self.run_time_option == 0:
            with open(os.path.join(self.data_dir, "ORIGINAL"    ), 'w') as f: pass
        elif self.run_time_option == 1:
            with open(os.path.join(self.data_dir, "QASM"), 'w') as f: pass
        else:
            with open(os.path.join(self.data_dir, "REAL"), 'w') as f: pass


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
                    Mrho[i] = self.Mrho_E[i]*( 1 + (-1 *self.ICrhoErrScale + 2*self.ICrhoErrScale*random.random()))

                    if Mrho[i] > 1: Mrho[i] = 1

                    Temp[i] = self.Temp_E[i]*( 1 + (-1 * self.ICtempErrScale + 2*self.ICtempErrScale*random.random()))

                    if Temp[i] > 1: Temp[i] = 1

            elif self.Shock_Flag == 1:
                if ((i == 0) or (i == ithroat - 1)):
                    Mrho[i] = self.Mrho_E[i]
                    Temp[i] = self.Temp_E[i]
        
                elif ( (i != 0) and (i != ithroat - 1) ):
                    Mrho[i] = self.Mrho_E[i]*( 1 + (-1 * self.ICrhoErrScale + 2*self.ICrhoErrScale*random.random()))

                    if Mrho[i] > 1: Mrho[i] = 1

                    Temp[i] = self.Temp_E[i]*( 1 + (-1 * self.ICtempErrScale + 2*self.ICtempErrScale*random.random()))

                    if Temp[i] > 1: Temp[i] = 1
                
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


    def IPrtn(self, b):
        #IPRTN partitions interval [Init, Fnl] into nn*NN subintervals
        hh = (b - self.a)/self.n     # width of a primary subinterval
        self.hbar = hh/self.N             # width of a sub-subinterval

        self.t = np.zeros((int(self.N), int(self.n)))        # create/initialize times array first (second)
        
        # slot specifies sub- (primary) subinterval
        # fill in t array-element values Note t(N,n) = Fnl
        for i in range(int(self.n)):
            for j in range(1, int(self.N) + 1):
                if i == 0: self.t[j - 1, i] = self.a + j * self.hbar    # fill times in first subinterval
                else: self.t[j - 1, i] = self.t[int(self.N) - 1,  i - 1] + j * self.hbar # fill remaining subintervals


    def IntegrateODE(self):
        #INTEGRATEODE numerically integrates ODE for 1D Navier-Stokes flow
        print("Running...\n")
        
        # Begin loop over the subintervals i result is approximate solution z(t).
        for i in prange(int(self.n)):
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

            # output initial simulation flow variables for subinterval i+1
            print('Code has completed subinterval ',i)

            if i != self.n: print('  Initial condition for next subinterval is:')
            elif i == self.n: print('  Final result for steady state U values are:')

            print(np.transpose(self.InitVal))

            print('Next subint start-time = ', (self.n**(self.k-1))*self.hbar*i)

            # if should log, then logging to the log file
            if self.log_QAmpEst:
                with open(self.QAmpEst_log_f, 'a') as f:
                    for item in self.data_vals: f.write(f"{item[0]}->{item[1]}\n")

                self.data_vals.clear()

            # write computational results to files for eventual plotting:
            print("Recording Results from Subinterval:", i)
            self.WriteResults()
            break

        print("\nDone...\n")


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


    def CalcFlux(self, U):
        #CALCFLUX evaluates flow fluxes at all grid-points for Nav-Stokes dynamics
        F = np.zeros((self.d, self.Tot_X_Pts))

        fac1 = (3 - self.Gamma)/2
        fac2 = (self.Gamma - 1)/self.Gamma
        fac3 = (self.Gamma - 1)/2

        # calculate fluxes loop over all grid-points 
        for i in range(self.Tot_X_Pts):
            F[0, i] = U[1, i]

            # introduce some useful definitions
            ratio1 = 1/U[0, i]
            ratio2 = (self.Gamma*U[1, i])/(U[0, i]*U[0, i])

            term1 = U[0, i]*U[2, i]
            term2 = U[1, i]*U[1, i]

            F[1, i] = ratio1*( fac2*term1 + fac1*term2 )
            F[2, i] = ratio2*( term1 - fac3*term2 )

        return F


    def CalcSource(self, U):
        #CALCSOURCE evaluates the source current in the Navier-Stokes dynamics
        J = np.zeros(self.Tot_Int_Pts)

        # calculate source current loop over interior grid-points with
        #       2 <= ll <= Tot_Int_Pts +1
        for i in range(1, (self.Tot_Int_Pts + 1)):    # ll indexes the grid-points                                       
            fac1 = ((self.Gamma - 1)/self.Gamma) * np.log(self.A[i+1]/self.A[i-1])
            fac2 = self.Gamma/2

            term1 = U[0, i] * U[2, i]
            term2 = U[1, i] * U[1, i]
                                    
            ratio1 = 1/U[0, i]

            J[i - 1] = fac1*ratio1*( term1 - fac2*term2 )


        return J


    def CalcFunc(self, F, J):
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

            ff_vals[0, i] = Fac1*( F[0, gridLabel + 1] - F[0, gridLabel - 1]) # HERE IS THE PROBLEM, THIS COLUMN IS NOT SET RIGHT
            ff_vals[1, i] = Fac1*( (F[1, gridLabel + 1] - F[1, gridLabel-1]) - J[i] )
            ff_vals[2, i] = Fac1*( F[2, gridLabel + 1] - F[2, gridLabel - 1] )


        return ff_vals


    def CalcfBvalsmSW(self, U, ff_vals):
        #CALCFBVALSMSW assigns ODE driver function boundary values - no shock-wave
        ff_Bvals = np.zeros((self.d,2))

        Fac3 = U[1, 0]/U[0, 0]

        # calculate ff_Bvals at nozzle entrance
        ff_Bvals[0, 0] = 0
        ff_Bvals[1, 0] = 2*ff_vals[1, 0] - ff_vals[1, 1]
        ff_Bvals[2, 0] = self.Gamma*Fac3*ff_Bvals[1, 0]

        # calculate use boundary conditions on flow variables to determine
        #   ff_Bvals at nozzle exit
        for k in range(3): ff_Bvals[k, 1] = 2*ff_vals[k, self.Tot_Int_Pts] - ff_vals[k, (self.Tot_Int_Pts-1)]

        return ff_Bvals


    def CalcfBvalspSW(self, U, ff_vals):
        #CALCFBVALSPSW assigns boundary values to ODE driver fcn - with shock-wave 
        ff_Bvals = np.zeros((self.d,2))

        Fac1 = self.Gamma/2
        Fac3 = U[1, 0]/U[0, 0]

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


    def Calc_dFdt(self, U, ff_vals, ff_Bvals):
        #CALC_DFDT evaluates first time derivative of flow fluxes at grid-points
        dFdt = np.zeros((self.d, self.Tot_X_Pts))

        fac4 = self.Gamma -1
        fac5 = self.Gamma/2
        fac6 = fac4*fac5
        fac7 = (3-self.Gamma)/2
        fac8 = fac4/self.Gamma

        # evaluate dFdt at boundaries:
        for ll in range(2):
            if ll == 0: # at nozzle entrance
                fac1 = U[1, 0]/U[0, 0]
                fac2 = U[2, 0]/U[0, 0]
                fac3 = fac1*fac2
                fac9 = (fac1)**(2)

                dFdt[0, ll] = ff_Bvals[2, ll]
                dFdt[1, ll] = fac7*fac1*(2*ff_Bvals[1, ll]-fac1*ff_Bvals[0, ll]) + fac8*ff_Bvals[2, ll]
                dFdt[2, ll] = self.Gamma*(fac1*ff_Bvals[2, ll]+fac2*ff_Bvals[1, ll]-fac3*ff_Bvals[0, ll])-fac6*fac9*(3*ff_Bvals[1, ll])-2*fac1*ff_Bvals[0, ll]

            elif ll == 1: # at nozzle exit
                fac1 = U[1, self.Tot_X_Pts - 1]/U[0, self.Tot_X_Pts - 1]
                fac2 = U[2, self.Tot_X_Pts - 1]/U[0, self.Tot_X_Pts - 1]
                fac3 = fac1*fac2
                fac9 = (fac1)**(2)

                dFdt[0, self.Tot_X_Pts- 1] = ff_Bvals[1, ll]
                dFdt[1, self.Tot_X_Pts- 1]= fac7*fac1*(2*ff_Bvals[1, ll]-fac1*ff_Bvals[0, ll]) + fac8*ff_Bvals[2, ll]
                dFdt[2, self.Tot_X_Pts- 1] = self.Gamma*(fac1*ff_Bvals[2, ll] + fac2*ff_Bvals[2, ll]- fac3*ff_Bvals[0, ll])-fac6*fac9*(3*ff_Bvals[1, ll]-2*fac1*ff_Bvals[0, ll])
            
            else:
                print('Unknown switch case label: ', ll)

        # evaluate dFdt at interior grid-points
        for ll in range(1, self.Tot_X_Pts):
            fac1 = U[1, ll]/U[0, ll]
            fac2 = U[2, ll]/U[0, ll]
            fac3 = fac1*fac2
            fac9 = (fac1)**(2)

            IPLabel = ll -1

            dFdt[0, ll] = ff_vals[1, IPLabel- 1]
            dFdt[1, ll] = fac7*fac1*(2*ff_vals[1, IPLabel- 1]-fac1*ff_vals[0, IPLabel- 1])+ fac8*ff_vals[2, IPLabel- 1]
            dFdt[2, ll] = self.Gamma*(fac1*ff_vals[2, IPLabel- 1] + fac2*ff_vals[1, IPLabel- 1] - fac3*ff_vals[1, IPLabel- 1])-fac6*fac9*(3*ff_vals[1, IPLabel- 1]-2*fac1*ff_vals[0, IPLabel- 1])

        return dFdt


    def Calc_dJdt(self, U, ff_vals):
        #CALC_DJDT evaluates first time derivative of source term at interior pts
        dJdt = np.zeros(self.Tot_Int_Pts)

        fac1 = (self.Gamma - 1)/self.Gamma
        fac2 = self.Gamma/2

        # evaluate dJdt by looping over interior grid-points 
        #           ( 2 <= ll < = self.Tot_Int_Pts + 1)
        for ll in range(self.Tot_Int_Pts + 1):
            IPLabel = ll - 1

            fac3 = U[1, ll]/U[0, ll]
            fac4 = np.log(self.A[ll+1]/self.A[ll-1])

            dJdt[IPLabel] = fac4*fac1*(ff_vals[2, IPLabel]-fac2*fac3*(2*ff_vals[1, IPLabel]-fac3*ff_vals[0, IPLabel] ) )

        return dJdt


    def Calc_dffdt(self, dFdt, dJdt):
        #CALC_DFFDT evaluates first time derivative of ODE driver function
        dffdt_vals = np.zeros((self.d, self.Tot_Int_Pts))
        fac1 = (-1)/(2*self.Del_x)

        # evaluate dffdt_vals at interior points (2 <= ll <= TotIPtp1)
        for ll in range(1, self.Tot_Int_Pts + 1):
            IPLabel = ll - 1

            dffdt_vals[0, IPLabel] = fac1*(dFdt[0, ll+1] - dFdt[0, ll-1])
            dffdt_vals[1, IPLabel] = fac1*(dFdt[1, ll+1] - dFdt[1, ll-1]-dJdt[IPLabel])
            dffdt_vals[2, IPLabel] = fac1*(dFdt[2, ll+1] - dFdt[2, ll-1])

        return dffdt_vals


    def CalcdfdtBvalspSW(self, U, ff_Bvals, dffdt_vals):
        #CALCDFDTBVALSPSW assigns dff/dt boundary values - shock-wave present
        dffdt_Bvals = np.zeros((self.d, 2))

        Fac1 = U[1, 0]/U[0, 0]
        Fac2 = (ff_Bvals[1, 0])**(2)/U[0, 0]

        # nozzle entrance calculations
        dffdt_Bvals[0, 0] = 0
        dffdt_Bvals[1, 0] = 2*dffdt_vals[1, 0] - dffdt_vals[1, 1]
        dffdt_Bvals[2, 0] = self.Gamma*( (Fac1*dffdt_Bvals[1, 0]) + Fac2 )

        # nozzle exit calculations
        Fac0 = self.Gamma/2
        Fac3 = U[1, self.Tot_X_Pts- 1]/U[0, self.Tot_X_Pts - 1]
        Fac4 = Fac3**(2)
        Fac5 = 1/U[0, self.Tot_X_Pts- 1]
        Fac6 = Fac3*Fac5
        Fac7 = Fac4*Fac5

        for k in range(self.d):
            dffdt_Bvals[k, 1] = 2*dffdt_vals[k, self.Tot_Int_Pts- 1]- dffdt_vals[k, (self.Tot_Int_Pts - 2)]

        dffdt_Bvals[2, 1] = Fac0*( Fac3*( 2*dffdt_Bvals[1, 1]- Fac3*dffdt_Bvals[0, 1] ) +2*Fac7*(ff_Bvals[0, 1])**(2) +2*Fac5*(ff_Bvals[1, 1])**(2)-4*Fac6*ff_Bvals[0, 1]*ff_Bvals[1, 1] )

        return dffdt_Bvals


    def CalcdfdtBvalsmSW(self, U, ff_Bvals, dffdt_vals):
        #CALCDFDTBVALSMSW assigns dff/dt boundary values - shock-wave absent
        dffdt_Bvals = np.zeros((self.d,2))

        Fac1 = U[1, 0]/U[0, 0]
        Fac2 = (ff_Bvals[1, 0])**(2)/U[0, 0]

        # nozzle entrance calculations
        dffdt_Bvals[0, 0] = 0
        dffdt_Bvals[1, 0] = 2*dffdt_vals[1, 0] - dffdt_vals[1, 1]
        dffdt_Bvals[2, 0] = self.Gamma*( (Fac1*dffdt_Bvals[1, 0]) + Fac2 )

        # nozzle exit calculations
        for k in range(3): dffdt_Bvals[k, 1] = 2*dffdt_vals[k, self.Tot_Int_Pts] - dffdt_vals[k, (self.Tot_Int_Pts - 1)]

        return dffdt_Bvals


    def Calc_d2Fdt2(self, U, dffdt_vals, dffdt_Bvals, ff_vals, ff_Bvals):
        #CALC_D2FDT2 evaluates flow flux second time derivative at all grid-points
        d2Fdt2 = np.zeros((self.d, self.Tot_X_Pts))

        fac4 = self.Gamma - 1
        fac5 = self.Gamma/2         
        fac6 = fac4*fac5           # (Gamma - 1)*(Gamma/2)
        fac7 = (3 - self.Gamma)/2
        fac8 = fac4/self.Gamma          # (Gamma - 1)/Gamma

        # evaluate d2Fdt2 at boundaries
        for ll in range(2):
            if ll == 0: # at nozzle ENTRANCE
                fac1 = U[1, ll]/U[0, ll]
                fac2 = U[2, ll]/U[0, ll]
                fac3 = fac1*fac2       # U(2,ll)U(3,ll)/(U(1,ll))**2
                fac9 = (fac1)**(2)      # (U(2,ll)/U(1,ll))**(2)
                fac10 = 1/U[0, ll]
                fac11 = fac9*fac10     # (U(2,ll))**(2)/(U(1,ll))**3
                fac12 = fac1*fac10     # U(2,ll)/(U(1,ll)**2
                fac13 = fac2*fac10     # U(3,ll)/(U(1,ll))**2
                fac14 = fac1*fac9      # (U(2,ll)/U(1,ll))**3
                fac15 = fac3*fac10     # (U(2,ll)U(3,ll))/(U(1,ll))**(3)
                fac16 = fac14*fac10    # (U(2,ll))**(3)/(U(1,ll))**(4)

                d2Fdt2[0, ll] = dffdt_Bvals[1, ll]
                d2Fdt2[1, ll] = fac7*( fac1*(2*dffdt_Bvals[1, ll]-fac1*dffdt_Bvals[0, ll])+2*fac10*(ff_Bvals[1, ll])**(2)+2*fac11*(ff_Bvals[0, ll])**(2)-4*fac12*ff_Bvals[0, ll]*ff_Bvals[1, ll])+fac8*dffdt_Bvals[2, ll]
                d2Fdt2[2, ll] = self.Gamma*( fac1*dffdt_Bvals[2, ll] + fac2*dffdt_Bvals[1, ll]-fac3*dffdt_Bvals[0, ll]+2*(fac10*ff_Bvals[1, ll]*ff_Bvals[2, ll]-fac12*ff_Bvals[0, ll]*ff_Bvals[2, ll]-fac13*ff_Bvals[0, ll]*ff_Bvals[1, ll]+fac15*(ff_Bvals[0, ll])**(2) ) ) - fac6*( 3*fac9*dffdt_Bvals[1, ll] -2*fac14*dffdt_Bvals[0, ll]+6*fac12*(ff_Bvals[1, ll])**(2) +6*fac16*(ff_Bvals[0, ll])**(2)-12*fac11*ff_Bvals[0, ll]*ff_Bvals[1, ll] ) 

            elif ll == 1:        # at nozzle EXIT
                fac1 = U[1, self.Tot_X_Pts- 1]/U[0, self.Tot_X_Pts- 1]
                fac2 = U[2, self.Tot_X_Pts- 1]/U[0, self.Tot_X_Pts- 1]
                fac3 = fac1*fac2       #U(2,TotXP)U(3,TotXPt)/(U(1,TotXPt))**2
                fac9 = (fac1)**(2)      # (U(2,TotXPt)/U(1,TotXPt))**(2)
                fac10 = 1/U[0, self.Tot_X_Pts- 1]
                fac11 = fac9*fac10     # (U(2,TotXPt))**(2)/(U(1,TotXPt))**3
                fac12 = fac1*fac10     # U(2,TotXPt)/(U(1,TotXPt)**2
                fac13 = fac2*fac10     # U(3,TotXPt)/(U(1,TotXPt))**2
                fac14 = fac1*fac9      # (U(2,TotXPt)/U(1,TotXPT))**3
                fac15 = fac3*fac10     #(U(2,ll)U(3,TotXPT))/(U(1,TotXPt))**(3)
                fac16 = fac14*fac10    # (U(2,TotXPt))**(3)/(U(1,TotXPt))**(4)

                d2Fdt2[0, self.Tot_X_Pts- 1] = dffdt_Bvals[1, ll]
                d2Fdt2[1, self.Tot_X_Pts- 1] = fac7*( fac1*(2*dffdt_Bvals[1, ll]-fac1*dffdt_Bvals[0, ll])+2*fac10*(ff_Bvals[1, ll])**(2)+2*fac11*(ff_Bvals[0, ll])**(2)-4*fac12*ff_Bvals[0, ll]*ff_Bvals[1, ll])+fac8*dffdt_Bvals[2, ll]
                d2Fdt2[2, self.Tot_X_Pts- 1] = self.Gamma*( fac1*dffdt_Bvals[2, ll] + fac2*dffdt_Bvals[1, ll]-fac3*dffdt_Bvals[0, ll]+2*(fac10*ff_Bvals[1, ll]*ff_Bvals[2, ll]-fac12*ff_Bvals[0, ll]*ff_Bvals[2, ll]-fac13*ff_Bvals[0, ll]*ff_Bvals[1, ll]+fac15*(ff_Bvals[0, ll])**(2) ) )-fac6*( 3*fac9*dffdt_Bvals[1, ll]-2*fac14*dffdt_Bvals[0, ll]+6*fac12*(ff_Bvals[1, ll])**(2)+6*fac16*(ff_Bvals[0, ll])**(2)-12*fac11*ff_Bvals[0, ll]*ff_Bvals[1, ll] )     
            
            else: print('Unknown switch case label: ', ll)

        # evaluate d2Fdt2 at interior points
        for ll in range(1, self.Tot_X_Pts - 1):
            IPLabel = ll - 1

            fac1 = U[1, ll]/U[0, ll]
            fac2 = U[2, ll]/U[0, ll]
            fac3 = fac1*fac2       # U(2,ll)U(3,ll)/(U(1,ll))**2
            fac9 = (fac1)**(2)      # (U(2,ll)/U(1,ll))**(2)
            fac10 = 1/U[0, ll]
            fac11 = fac9*fac10     # (U(2,ll))**(2)/(U(1,ll))**3
            fac12 = fac1*fac10     # U(2,ll)/(U(1,ll)**2
            fac13 = fac2*fac10     # U(3,ll)/(U(1,ll))**2
            fac14 = fac1*fac9      # (U(2,ll)/U(1,ll))**3
            fac15 = fac3*fac10     # (U(2,ll)U(3,ll))/(U(1,ll))**(3)
            fac16 = fac14*fac10    # (U(2,ll))**(3)/(U(1,ll))**(4)

            d2Fdt2[0, ll] = dffdt_vals[1, IPLabel- 1]
            d2Fdt2[1, ll] = fac7*( fac1*(2*dffdt_vals[1, IPLabel- 1]-fac1*dffdt_vals[0, IPLabel- 1])+2*fac10*(ff_vals[1, IPLabel- 1])**(2)+2*fac11*(ff_vals[0, IPLabel- 1])**(2)-4*fac12*ff_vals[0, IPLabel- 1]*ff_vals[1, IPLabel- 1])+fac8*dffdt_vals[2, IPLabel- 1]
            d2Fdt2[2, ll] = self.Gamma*( fac1*dffdt_vals[2, IPLabel- 1] + fac2*dffdt_vals[1, IPLabel- 1]-fac3*dffdt_vals[0, IPLabel- 1]+2*(fac10*ff_vals[1, IPLabel- 1]*ff_vals[2, IPLabel- 1]-fac12*ff_vals[0, IPLabel- 1]*ff_vals[2, IPLabel- 1]-fac13*ff_vals[0, IPLabel- 1]*ff_vals[1, IPLabel- 1]+fac15*(ff_vals[0, IPLabel- 1])**(2) ) )-fac6*( 3*fac9*dffdt_vals[1, IPLabel- 1]-2*fac14*dffdt_vals[0, IPLabel- 1]+6*fac12*(ff_vals[1, IPLabel- 1])**(2)+6*fac16*(ff_vals[0, IPLabel- 1])**(2)-12*fac11*ff_vals[0, IPLabel- 1]*ff_vals[1, IPLabel- 1] )

        return d2Fdt2


    def Calc_d2Jdt2(self, U, dffdt_vals, ff_vals):
        #CALC_D2JDT2 evaluates 2nd time derivative of source term at interior pts
        d2Jdt2 = np.zeros((self.Tot_Int_Pts))

        fac1 = (self.Gamma - 1)/self.Gamma
        fac2 = self.Gamma/2

        # evaluate d2Jdt2 by looping over interior grid-points
        #           ( 2 <= ll <= TotIPtp1 )
        for ll in range(1, self.Tot_Int_Pts + 1):
            IPLabel = ll -1

            fac3 = U[1, ll]/U[0, ll]
            fac4 = np.log(self.A[ll+1]/self.A[ll-1])
            fac5 = (fac1)**(2)          # (U(2,ll)/U(1,ll))**(2)
            fac6 = 1/U[0, ll]           # 1/U(1,ll)
            fac7 = (fac6)**(2)          # 1/(U(1,ll))**(2)
            fac8 = fac6*fac7           # 1/(U(1,ll))**(3)
            fac9 = fac3*fac6           # U(2,ll)/(U(1,ll))**(2)

            d2Jdt2[IPLabel] = fac4*fac1*(dffdt_vals[2, IPLabel] -fac2*( 2*fac3*dffdt_vals[1, IPLabel] -fac5*dffdt_vals[0, IPLabel]+2*fac6*(ff_vals[1, IPLabel])**(2)+2*fac8*(U[1, ll]*ff_vals[0, IPLabel])**(2) - 4*fac9*ff_vals[0, IPLabel]*ff_vals[1, IPLabel] ) )

        return d2Jdt2


    def Calc_d2ffdt2(self, d2Fdt2, d2Jdt2):
        #CALC_D2FFDT2 evaluates second time derivative of ODE driver function
        d2ffdt2_vals = np.zeros((self.d, self.Tot_Int_Pts))

        fac1 = (-1)/(2*self.Del_x)

        # evaluate d2ffdt2_vals at interior grid-points
        #           ( 2 <= ll <= TotIPtp1 )
        for ll in range(1, self.Tot_Int_Pts + 1):
            IPLabel = ll -1  

            d2ffdt2_vals[0, IPLabel] = fac1*(d2Fdt2[0, ll+1] - d2Fdt2[0, ll-1])
            d2ffdt2_vals[1, IPLabel] = fac1*(d2Fdt2[1, ll+1] - d2Fdt2[1, ll-1]-d2Jdt2[IPLabel])                
            d2ffdt2_vals[2, IPLabel] = fac1*(d2Fdt2[2, ll+1] - d2Fdt2[2, ll-1])    

        return d2ffdt2_vals


    def BldTMat(self, f):
        #BLDTMAT constructs Taylor polynomial coefficients 
        rmax = self.r + 1       # degree of Taylor polynomials
        rmaxp1 = rmax + 1  # number of terms in Taylor polynomials

        # for given subsubinterval mat will return the Taylor polynomial 
        #  coefficients for each component of approximation to ODE driver 
        #  function and interior grid-point
        mat = np.zeros((self.d, rmaxp1, self.Tot_Int_Pts))

        # evaluate Taylor polynomial coefficients by looping around the interior
        #       grid-points, then the d components, then the rmax coefficients
        for ll in range(self.Tot_Int_Pts):
            for k in range(self.d):
                for m in range(rmax): mat[k, m, ll] = f[k, (rmax - m - 1), ll]/math.factorial(rmax - m - 1)

                # assign value of polynomial at start of subsubinterval. since InVal
                #     is defined at all grid-points, grid-point label is ll + 1 for
                #     InVal
                mat[k, rmaxp1- 1, ll] = self.InitVal[k, (ll+1)]

        return mat


    def NextInCond(self, U, mmat):
        #NEXTINCOND determines initial condition for next subsubinterval 
        rmax = self.r + 1       # degree of Taylor Polynomials
        rmaxp1 = rmax + 1  # number of terms in Taylor polynomials

        TylrPoly = np.zeros(rmaxp1) # array to store coefficients for a Taylor
            #   polynomial
            
        NextInVal = np.zeros((self.d, self.Tot_X_Pts))     # array to store initial conditions 
                    #   for next subsubinterval 

        # calculate initial condition for next subsubinterval for the
        #       interior grid-points    ( 2 <= ll <= TotXPtm1 )
        for ll in range(self.Tot_Int_Pts):
            GPLabel = ll + 1

            for pp in range(self.d):
                for mm in range(rmax): TylrPoly[mm] = mmat[pp, mm, ll]

                TylrPoly[rmaxp1- 1] = U[pp, GPLabel]
                NextInVal[pp, GPLabel] = np.polyval(TylrPoly, self.hbar)

        # use flow boundary conditions to determine new initial conditions for
        #       boundary grid-points
        if self.Shock_Flag == 0: U_Bvals = self.CalcBCmSW(U)
        elif self.Shock_Flag == 1: U_Bvals = self.CalcBCpSW(U)
        else: print(' Unknown Shock_Flag value: ', self.Shock_Flag)

        # store boundary values in NextInVal
        for p in range(self.d):
            NextInVal[p, 0] = U_Bvals[p, 0]      # values at nozzle entrance
            NextInVal[p, self.Tot_X_Pts - 1] = U_Bvals[p, 1]  # values at nozzle exit

        return NextInVal


    def CalcBCpSW(self, U):
        #CALCBCPSW evaluates flow variables at nozzle boundaries-shock-wave present
        U_Bvals = np.zeros((self.d, 2))

        fac1 = 1/(self.Gamma -1)
        fac2 = self.Gamma/2

        # evaluate U at nozzle entrance:
        U_Bvals[0, 0] = self.A[0]
        U_Bvals[1, 0] = 2 * U[1, 1] - U[1, 2]

        fac3 = (U_Bvals[1, 0]/U_Bvals[0, 0])**(2)

        U_Bvals[2, 0] = U_Bvals[0, 0]*(fac1 + fac2*fac3)

        # evaluate flow at nozzle exit:
        for k in range(self.d - 1):
            U_Bvals[k, 1] = 2 * U[k, self.Tot_X_Pts - 2] - U[k, self.Tot_X_Pts - 3]

        U_Bvals[2, 1] = self.Exit_Pressure*self.A[self.Tot_X_Pts- 1]*fac1 +( fac2*( U_Bvals[1, 1])**(2) )/U_Bvals[0, 1]

        return U_Bvals


    def CalcBCmSW(self, U):
        #CALCBCMSW evaluates flow variables at nozzle boundaries-shock-wave absent
        U_Bvals = np.zeros(self.d,2)

        fac1 = 1/(self.Gamma -1)
        fac2 = self.Gamma/2

        # evaluate U at nozzle entrance:
        U_Bvals[0, 0] = self.A[0]
        U_Bvals[1, 0] = 2*U[1, 1] - U[1, 2]

        fac3 = (U_Bvals[1, 0]/U_Bvals[0, 0])**(2)

        U_Bvals[2, 0] = U_Bvals[0, 0]*(fac1 + fac2*fac3)

        # evaluate flow at nozzle exit:
        for k in range(self.d): U_Bvals[k, 1] = 2 * U[k, self.Tot_X_Pts -2] - U[k, self.Tot_X_Pts -3]

        return U_Bvals


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
        for j in range(self.N):
            #start = time.time()
            print('In subinterval i =', i, '   starting subsubinterval j =', j)

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
            for k in range(self.d):
                for ll in range(self.Tot_Int_Pts):
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

        # Eq. (28) in Kacewicz requires Gint above to be divided by N ( = ml)
        #  other division by N included in calculation of mean by MeanOrc
        return Gint/self.N


    def fOrc(self, t, Start, TCoeffs):
        #FORC evaluates ODE driver function f at l**{s}_{i}(t) at interior grd-pts
        # initialize parameter and array
        delt = t - Start            # elapsed time from start of subsubinterval j
        lt = np.zeros((self.d, self.Tot_Int_Pts))   # stores l**{s}_{i}(t) at each interior 
                            #     grid-point

        Poly = np.zeros(int(self.r + 2))  # initialize to zero array storing Taylor 
                        #    polynomial coefficients for l(t(j-1), i) 
                        #    for given component and interior grid-point
                        
        U = np.zeros((self.d, self.Tot_X_Pts))  # array to store primary flow variables

        # evaluate l**{s}_{i}(t) at each interior grid-point, one component at time
        for ll in range(self.Tot_Int_Pts):    
            for m in range(self.d):
                # load Taylor polynomial coefficients for m-th component at
                #   interior grid-point ll in array Poly
                for column in range(self.r + 2): Poly[column] = TCoeffs[m, column, ll]

                # store value of m-th Taylor polynomial at elapsed time delt and 
                #       at interior grid-point ll
                lt[m, ll] = np.polyval(Poly, delt) 


        # assign U at each interior grid-point
        for ll in range(1, self.Tot_X_Pts - 1):
            for m in range(self.d): U[m, ll] = lt[m, ll - 1]


        # assign U at boundary points using flow boundary conditions
        if self.Shock_Flag == 0: U_Bvals = self.CalcBCmSW(U)
        elif self.Shock_Flag == 1: U_Bvals = self.CalcBCpSW(U)
        else: print('Unknown Shock_Flag value: ', self.Shock_Flag)

        for m in range(self.d):
            U[m, 0] = U_Bvals[m, 0]
            U[m, self.Tot_X_Pts - 1] = U_Bvals[m, 1]

        #  evaluate f0 which stores ODE driver function values at all interior
        #   grid-points at start of subsubinterval j. f0 is d x Tot_Int_Pts array
        return self.CalcFunc(self.CalcFlux(U), self.CalcSource(U))


    def QAmpEst(self, omega):
        # runs the function that was parsed in InitCalcParms
        result = self.algo_func(omega)
        
        self.data_vals.append([omega, result])

        return result

    
    def original_algo(self, omega):
        # run in the c++ file
        return self.c_lib(ctypes.c_int(self.N), ctypes.c_double(self.delta1), ctypes.c_double(omega))
    

    def simulator_algo(self, omega):
        # run on the qasm simulator
        A = BernoulliA(omega)
        Q = BernoulliQ(omega)

        problem = EstimationProblem(
            state_preparation=A,  # A operator
            grover_operator=Q,  # Q operator
            objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
        )

        return (np.sin(np.pi * self.ae.estimate(problem).estimation))**2 #(np.sin(np.pi * self.ae.estimate(problem).mle))**2


    def real_algo(self, omega):
        A = BernoulliA(omega)
        Q = BernoulliQ(omega)

        problem = EstimationProblem(
            state_preparation=A,  # A operator
            grover_operator=Q,  # Q operator
            objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
        )

        return (np.sin(np.pi * self.ae.estimate(problem).mle))**2


    def Calc_FlowVarResults(self):
        #CALC_FLOWVARRESULTS finds physical flow variables from simulation results
        # define useful factors
        fac1 = self.Gamma - 1
        fac2 = self.Gamma/2

        # loop over all grid-points
        for i in range(self.Tot_X_Pts):
            self.Mrho_D[i] = self.InitVal[0, i] / self.A[i]
            self.Vel_D[i] = self.InitVal[1, i] / self.InitVal[0, i]
            self.Temp_D[i] = fac1*((self.InitVal[2, i] / self.InitVal[0, i]) - fac2 * self.Vel_D[i] * self.Vel_D[i])
            self.Press_D[i] = self.Mrho_D[i] * self.Temp_D[i]
            self.Mach_D[i] = self.Vel_D[i] / np.sqrt(self.Temp_D[i])


    def WriteResults(self):       
        #WRITERESULTS writes results of Navier-Stokes solution to files.
        # calculate relative error in primary flow variables
        Rel_MachErr = np.divide(np.abs(self.Mach_D - self.Mach_E), self.Mach_E)
        Rel_MrhoErr = np.divide(np.abs(self.Mrho_D - self.Mrho_E), self.Mrho_E)
        Rel_PressErr = np.divide(np.abs(self.Press_D - self.Press_E), self.Press_E)
        Rel_TempErr = np.divide(np.abs(self.Temp_D - self.Temp_E), self.Temp_E)
        Rel_VelErr = np.divide(np.abs(self.Vel_D - self.Vel_E), self.Vel_E)

        # calculate mean mass flow rate at final time define array to store
        MeanU2 = np.mean(self.U2[: , int(self.n)])

        # calculate mean and standard deviation of relative errors and store
        MeanRelTempErr = np.mean(Rel_TempErr)
        MeanRelMachErr = np.mean(Rel_MachErr)
        MeanRelMrhoErr = np.mean(Rel_MrhoErr)
        MeanRelPressErr = np.mean(Rel_PressErr)
        MeanRelVelErr = np.mean(Rel_VelErr)

        SDevRelTempErr = np.std(Rel_TempErr)
        SDevRelMachErr = np.std(Rel_MachErr)
        SDevRelMrhoErr = np.std(Rel_MrhoErr)
        SDevRelPressErr = np.std(Rel_PressErr)
        SDevRelVelErr = np.std(Rel_VelErr)

        files = [
                "self.U2",
                "MeanU2",
                "MeanRelTempErr",
                "MeanRelMachErr",
                "MeanRelMrhoErr",
                "MeanRelPressErr",
                "MeanRelVelErr",
                "SDevRelTempErr",
                "SDevRelMachErr",
                "SDevRelMrhoErr",
                "SDevRelPressErr",
                "SDevRelVelErr",
                "self.Mach_D",
                "self.Mach_E",
                "self.Mrho_D",
                "self.Mrho_E",
                "self.Press_D",
                "self.Press_E",
                "self.Temp_D",
                "self.Temp_E",
                "self.ff0_throat",
                "self.ff1_throat",
                "self.ff2_throat",
                "Rel_MachErr",
                "Rel_MrhoErr",
                "Rel_PressErr",
                "Rel_TempErr",
                "Rel_VelErr"
                ]


        for file in files:
            val = eval(file)
            file = os.path.join(self.data_dir, file.replace("self.", ''))
            if isinstance(val, np.ndarray): np.savetxt(file, val, "%8.6f")
            else:
                with open(file, 'w') as f: f.write(str(val))


@jit
def run():
    inst = ns_q()

run()

