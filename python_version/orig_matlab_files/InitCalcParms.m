function [Del_x, A, Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, ...
          InitVal, n, N, t, hbar, Tot_Int_Pts, U2, delta1, k, ithroat, ...
          ff0_throat, ff1_throat, ff2_throat] = ...
                InitCalcParms(x_min, x_max, Tot_X_Pts, Gamma, ...
                              Exit_Pressure, Shock_Flag, In_Mass_Flow, ...
                              d, ICMFlowErrScale, ICrhoErrScale, ...
                              ICtempErrScale, in_n, delta, err1, ...
                              Tot_TSteps, a)
%INITCALCPARMS initializes inputs to quantum 1D Navier-Stokes algorithm.
%   InitCalcParms initializes all the inputs needed by quantum algorithm 
%       to solve 1D Navier-Stokes equations.
%
%  INPUTS:
%       x_min = initial spatial grid-point/nozzle entrance
%       x_max = final spatial grid-point/nozzle exit
%       Tot_X_Pts = total number of spatial grid-points used
%       Gamma = ratio of specific heats c_p/c_v
%       Exit_Pressure = pressure at nozzle exit
%       Shock_Flag = 0 (1) if shock-wave absent (present)
%       In_Mass_Flow = initial value for mass flow through nozzle
%       d = number of components of ODE solution z(t)
%       ICMFlowErrScale = scale of random shift to initial mass flow
%       ICrhoErrScale = scale of random shift to initial mass density
%       ICtempErrScale = scale of random shift to initial temperature
%       in_n = initial (guess of) number of time subintervals in [a,b]
%       delta = probability that approximate solution z(t) violates its
%                error upper bound (Eq. (13) in Kacewicz quantum ODE paper)
%       err1 = upper bound on error of approximate value found for function
%               mean using quantum amplitude estimation algorithm
%       Tot_TSteps = total number of time-steps taken by solver
%       a = initial time
%
%  OUTPUTS:
%       Del_x = spatial separation of adjacent grid-points
%       A = 1 x Tot_X_Pts array storing nozzle area at each grid-point
%       Mach_E = 1 x Tot_X_Pts array storing exact result for steady-state
%                   Mach number at all grid-points
%       Mrho_E = 1 x Tot_X_Pts array storing exact result for steady-state
%                   mass density at all grid-points
%       Press_E = 1 x Tot_X_Pts array storing exact result for steady-state
%                   pressure at all grid-points
%       Temp_E = 1 x Tot_X_Pts array storing exact result for steady-state
%                   temperature at all grid-points
%       Vel_E = 1 x Tot_X_Pts array storing exact result for steady-state
%                   velocity at all grid-points
%       SW_JumpP = ratio of pressure behind SW/pressure ahead of SW
%       SW_JumpM = ratio of Mach number behind SW/Mach number ahead of SW
%       InitVal = d x Tot_X_Pts array storing initial values of
%                   computational flow variables U(d,Tot_X_Pts)
%       n = number of time subintervals in [a,b] used in simulation;
%               insures hbar is less than CFL time-step upper bound.
%       N = number of time sub-subintervals used in simulation
%       t = N x n array storing partition times of subsubintervals j inside
%               subinterval i. Specifically, t(j,i) is largest time in
%               sub-subinterval j in subinterval i. Note that smallest time
%               in sub-subinterval j is t(j-1,i) and the initial time a is
%               NOT stored in t.
%       hbar = width of a time sub-subinterval
%       Tot_Int_Pts = number of interior grid-points
%       U2 = Tot_X_Pts x (n+1) array which stores the calculated mass flow
%               rate at beginning of each time subinterval as well as the
%               mass flow rate at the final time of the ODE integration
%       delta1 = probability that approximate value of function mean
%                   violates its upper bound
%       k = number of recursion levels used
%       ithroat = grid-point index for nozzle throat
%       ff0_throat = driver function value at throat at end of each
%                       time subinterval
%       ff1_throat = value of first derivative of driver function at 
%                       throat at end of each time subinterval
%       ff2_throat = value of second derivative of driver function at 
%                       throat at end of each time subinterval
%
%
%  Support functions: Set_XGrid, Calc_Noz_Area, Calc_ExactResultsmSW,
%                       Calc_ExactResultspSW, SetInCond, InitParms,
%                           IPrtn

% set up x-grid for problem

  [x, Del_x] = Set_XGrid(x_min, x_max, Tot_X_Pts);

  %x
  %Del_x
  
% calculate nozzle area at each grid-point

  A = Calc_Noz_Area(x);

  %A
  
% calculate exact steady-state solution of primary flow variables: Mach
%   number, mass density, pressure, and temperature

  if Shock_Flag == 0
      [Mach_E, Mrho_E, Press_E, Temp_E, Vel_E] = ...
          Calc_ExactResultsmSW(Gamma, Tot_X_Pts, A);
  elseif Shock_Flag == 1
      [Mach_E, Mrho_E, Press_E, Temp_E, Vel_E, SW_JumpP, SW_JumpM] ...
        = Calc_ExactResultspSW(Gamma, Tot_X_Pts, A, Exit_Pressure, Del_x);
  else
      display(['Unknown Shock_Flag value: ' int2str(Shock_Flag)])
  end


% Initialize initial condition InitVal of ODE solution. Delta_t is maximum
%   time-step satisfying Courant-Friedrichs-Levy (CFL) numerical stability
%   condition. InitVal is d x Tot_X_Pts array storing initial 
%   computational flow variables U(d,Tot_X_Pts)

%Shock_Flag
%In_Mass_Flow
%Gamma
%x
%Del_x
%A
%d
%Mrho_E
%Temp_E
%ICMFlowErrScale
%ICrhoErrScale
%ICtempErrScale

  [InitVal, Delta_t, In_Mass_Flow_Noisy] = ...
     SetInCond(Shock_Flag, In_Mass_Flow, Gamma, x, Del_x, A, d, Mrho_E,...
               Temp_E, ICMFlowErrScale, ICrhoErrScale, ICtempErrScale);

%InitVal
% calculate final time b

  b = Tot_TSteps*Delta_t;
  
% initialize parameters for quantum ODE algorithm

  [N, delta1, final_n, k] = InitParms(in_n, delta, err1, Tot_TSteps);
  
  n = final_n;
  
% Partition time interval [a, b] into n^(k) sub-subintervals by introducing
%   N - 1 intermediate times t(j,i), where t(j,i) is largest time in
%   sub-subinterval j in subinterval i. The smallest time in
%   sub-subinterval j is t(j-1,i). Note that the initial time t = a is NOT
%   stored in N x n array t

  [t, hbar] = IPrtn(a, b, n, N);
  
% Taylor polynomials only needed at interior grid-points so define

  Tot_Int_Pts = Tot_X_Pts - 2;
  
% initialize array U2 to store calculated mass flow rate at beginning of
%   each subinterval i and also at time at end of ODE integration.
%       U2 = Tot_X_Pts x (n+1) array

  U2 = zeros(Tot_X_Pts, n+1);
  
% use initial condition to assign U2 at start of subinterval i = 1

  for gridpt = 1:Tot_X_Pts
      U2(gridpt, 1) = InitVal(2, gridpt);
  end
  
% define ithroat as grid-point index of nozzle throat

  ithroat = (Tot_X_Pts + 1)/2;
  
% initialize arrays to store value of driver function f(z(t)) and its first
%   r time derivatives at nozzle throat at end of each subinterval i

  ff0_throat = zeros(d,n);
  ff1_throat = zeros(d,n);
  ff2_throat = zeros(d,n);

  
% output initial condition for runtime inspection

  %usr_message = 'Calculation paused - press any key to continue: ';

  InitVal'
  %disp('Above are initial values of computational flow values U')
  
  %disp(usr_message)
  %pause
  
% output following parameters for runtime inspection

  Delta_t
  hbar
  b
  n
  k
  In_Mass_Flow_Noisy
  
  if Shock_Flag == 1
      SW_JumpP
      SW_JumpM
  end


  %disp(usr_message)
  %pause
  

% start timer for the calculation in driver once this function returns

end

