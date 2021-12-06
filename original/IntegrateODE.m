function [U2, Mach_D, Mrho_D, Press_D, Temp_D, Vel_D, Rel_MachErr, ...
          Rel_MrhoErr, Rel_PressErr, Rel_TempErr, Rel_VelErr, ...
          AvRelTempErr, AvPlusSDevRelTempErr, AvMinusSDevRelTempErr,...
          AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr,...
          AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr,...
          AvRelPressErr, AvPlusSDevRelPressErr, AvMinusSDevRelPressErr,...
          AvU2, ...
          ff0_throat, ff1_throat, ff2_throat] ...
    = IntegrateODE(d, n, N, hbar, r, Del_x, Gamma, Tot_Int_Pts, k, ...
                   Tot_X_Pts, Shock_Flag, Exit_Pressure, ithroat, a,...
                   delta1, rho, InitVal, A, t, U2_in, ff0_throat_in,...
                   ff1_throat_in, ff2_throat_in, Mach_E, Mrho_E, ...
                   Press_E, Temp_E, Vel_E, In_Mass_Flow)
%INTEGRATEODE numerically integrates ODE for 1D Navier-Stokes flow
%   IntegrateODE numerically integrates the set of ordinary differential
%       equations for 1D Navier-Stokes flow using Kacewicz' quantum
%       algorithm: B. Kacewicz, J. Complexity vol. 22 (2006) 676-690.
%
%   INPUTS:
%       d = number of components of ODE solution z(t)
%       n = number of subintervals used in time partition
%       N = number of time sub-subintervals used in integration
%       hbar = width of a time sub-subinterval
%       r = number of derivatives allowed for ODE driver function f[z]
%       Del_x = spatial separation of adjacent grid-points
%       Gamma = ratio of specific heats c_p/c_v
%       Tot_Int_Pts = total number of interior grid-points
%       k = number of recursion levels used
%       Tot_X_Pts = total number of spatial grid-points
%       Shock_Flag = 0 (1) when shock-wave is absent (present)
%       Exit_Pressure = dimensionless pressure at nozzle exit
%       ithroat = grid-point index of nozzle throat
%       a = initial time for integration
%       delta1 = probability that approximate value of function mean
%                   violates its upper bound
%       rho = basic Holder class parameter; tentatively set equal to 1
%       InitVal = d x Tot_X_Pts array storing initial values of
%                   computational flow variables for time subinterval -
%                   assign value for start of subinterval i+1 at end of
%                   subinterval i integration
%       A = 1 x Tot_X_Pts array storing dimensionless nozzle area at each
%               grid-point
%       t = N x n array storing partition times of sub-subintervals j
%              inside subinterval i. Specifically, t(j,i) is largest time
%              in sub-subinterval j in subinterval i. Note that smallest
%              time in sub-subinterval j is t(j-1,i) and initial time a is
%              NOT stored in t.
%       U2_in = Tot_X_Pts x (n+1) array storing initial mass flow rate at
%                   all grid-points in column 1 - all other columns are
%                   zero. This is immediately assigned to U2 and
%                   integration of ODE fills remaining columns
%       ff0_throat_in = d x n array storing driver function value at 
%                           throat; value passed is the all zeros array
%       ff1_throat_in = d x n array storing first derivative of driver
%                           function at throat; value pass is zero array
%       ff2_throat_in = d x n array storing second derivative of driver
%                           function at throat; value pass is zero array
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
%       In_Mass_Flow = initial mass flow rate through nozzle
%
%   OUTPUTS:
%       U2 = Tot_X_Pts x (n+1) array storing calculated mass flow rate at
%               beginning of each subinterval, and at the final time of ODE
%               integration
%       Mach_D = 1 x Tot_X_Pts array storing calculated Mach Number at each
%                   grid-point at end of ODE integration
%       Mrho_D = 1 x Tot_X_Pts array storing calculated dimensionless mass
%                   density at each grid-point at end of ODE integration
%       Press_D = 1 x Tot_X_Pts array storing calculated dimensionless 
%                   pressure at each grid-point at end of ODE integration
%       Temp_D = 1 x Tot_X_Pts array storing calculated dimensionless 
%                   temperature at each grid-point at end of ODE
%                   integration
%       Vel_D = 1 x Tot_X_Pts array storing calculated dimensionless 
%                   velocity at each grid-point at end of ODE
%                   integration
%       Rel_MachErr = 1 x Tot_X_Pts array storing relative error in Mach
%                       number at each grid-point at end of ODE integration
%       Rel_MrhoErr = 1 x Tot_X_Pts array storing relative error in
%                       dimensionless mass density at each grid-point at 
%                       end of ODE integration
%       Rel_PressErr = 1 x Tot_X_Pts array storing relative error in
%                       dimensionless pressure at each grid-point at 
%                       end of ODE integration
%       Rel_TempErr = 1 x Tot_X_Pts array storing relative error in
%                       dimensionless temperature at each grid-point at 
%                       end of ODE integration
%       Rel_VelErr = 1 x Tot_X_Pts array storing relative error in
%                       dimensionless velocity at each grid-point at 
%                       end of ODE integration
%       AvRelTempErr = 1 x Tot_X_Pts array storing average relative
%                       temperature error, same value to each grid-point
%       AvPlusSDevRelTempErr = 1 x Tot_X_Pts array storing (average plus
%                       standard deviation) relative temperature error, 
%                       same value to each grid-point
%       AvMinusSDevRelTempErr = 1 x Tot_X_Pts array storing (average minus
%                       standard deviation) relative temperature error, 
%                       same value to each grid-point
%       AvRelMachErr = 1 x Tot_X_Pts array storing average relative Mach
%                       number error, same value to each grid-point
%       AvPlusSDevRelMachErr = 1 x Tot_X_Pts array storing (average plus
%                       standard deviation) relative Mach number error, 
%                       same value to each grid-point
%       AvMinusSDevRelMachErr = 1 x Tot_X_Pts array storing (average minus
%                       standard deviation) relative Mach number error, 
%                       same value to each grid-point
%       AvRelMrhoErr = 1 x Tot_X_Pts array storing average relative mass
%                       density error, same value to each grid-point
%       AvPlusSDevRelMrhoErr = 1 x Tot_X_Pts array storing (average plus
%                       standard deviation) relative mass density error, 
%                       same value to each grid-point
%       AvMinusSDevRelMrhoErr = 1 x Tot_X_Pts array storing (average minus
%                       standard deviation) relative mass density error, 
%                       same value to each grid-point
%       AvRelPressErr = 1 x Tot_X_Pts array storing average relative 
%                       pressure error, same value to each grid-point
%       AvPlusSDevRelPressErr = 1 x Tot_X_Pts array storing (average plus
%                       standard deviation) relative pressure error, 
%                       same value to each grid-point
%       AvMinusSDevRelPressErr = 1 x Tot_X_Pts array storing (average minus
%                       standard deviation) relative pressure error, 
%                       same value to each grid-point
%       AvU2 = 1 x Tot_X_Pts array storing average mass flow rate (over 
%                 nozzle) at final time of ODE integration at each grid-pt
%       ff0_throat = d x n array storing driver function value at throat at
%                       end of each subinterval
%       ff1_throat = d x n array storing first derivative of driver 
%                      function value at throat at end of each subinterval
%       ff2_throat = d x n array storing second derivative of driver 
%                      function value at throat at end of each subinterval
%
%   Support functions: BldTPoly, IntegrateGij, CalcBCmSW, CalcBCpSW,
%                       Calc_FlowVarResults

% initialize arrays ff0_throat, ff1_throat, ff2_throat to values passed
%   through function input arguments 

 ff0_throat = ff0_throat_in;
 ff1_throat = ff1_throat_in;
 ff2_throat = ff2_throat_in;
 
% similarly, initialize U2 to U2_in

 U2 = U2_in;

% Begin loop over the subintervals i; result is approximate solution z(t).

 for i = 1:n
    %build Taylor polynomials l^{s}_{i}(t)for subinterval i at all
    %   interior grid-points; store polynomial coefficients in 
    %       StoreLz(d,r+2,Tot_Int_Pts,N)
    %   ff_throat is a d x (r+1) array storing the values of residual and
    %       its first r time derivatives at the nozzle throat at the end 
    %       of subinterval i
    
    [StoreLz, ff_throat] = BldTPoly(d,n,N,hbar,r,InitVal,Del_x,...
                                        Gamma,Tot_Int_Pts,...
                                            Tot_X_Pts, A, Shock_Flag,...
                                                Exit_Pressure, ithroat);
    
    % store values of ff_throat for subinterval i for easy of writing to
    %   files
    
    ff0_throat( :, i) = abs(ff_throat( :, 1));
    ff1_throat( :, i) = abs(ff_throat( :, 2));
    ff2_throat( :, i) = abs(ff_throat( :, 3));
                            
    % store N intermediate times for subinterval i
    
    StoreTimes4i = t(:,i);  
    
    % define Start time for subinterval i
    
    if i == 1
        Start = a;
    else
        Start = t(N,i-1);
    end
    
    %  gInt = d x Tot_Int_Pts array storing integral of each component of
    %               g_ij over subinterval i for each interior grid-point
    %
    
    [gInt] = IntegrateGij(StoreLz, StoreTimes4i, Start, ...
                                      d, r, N, delta1, hbar, rho, ...
                                        Tot_Int_Pts, A, Gamma, Del_x, ...
                                          Shock_Flag, i, Exit_Pressure);
         
    % add gInt for subint i to InitVal to get InitVal for 
    %   subinterval i + 1 at each interior grid-point
    
    for ll = 2:(Tot_Int_Pts+1)
        IP_Label = ll - 1;
        
        for m = 1:d
            InitVal(m,ll) = InitVal(m,ll) + gInt(m,IP_Label);
        end
    end
    
    % use flow boundary conditions to determine InitVal for
    %   subinterval i + 1 at two boundary points 1 and Tot_X_Pts. Store
    %   them in InitVal_BVals
    
    if Shock_Flag == 0
        InitVal_BVals = CalcBCmSW(InitVal,A,Gamma,d,Tot_X_Pts);
    elseif Shock_Flag == 1
        InitVal_BVals = CalcBCpSW(InitVal,A,Gamma,d,Tot_X_Pts,...
                                    Exit_Pressure);
    else
        disp(['Unknown Shock_Flag value: ' int2str(Shock_Flag)])
    end
    
    % transfer boundary values from InitVal_Bvals to InitVal. Resulting 
    %       InitVal then contains initial values of primary flow 
    %       variable U for subinterval i + 1 at each grid-point.
    
    for m = 1:d
        InitVal(m,1) = InitVal_BVals(m,1);
        InitVal(m,Tot_X_Pts) = InitVal_BVals(m,2);
    end
    
    % store U2 values at start of subinterval i + 1 at all gridpoints
     
    for gridpt = 1:Tot_X_Pts
        U2(gridpt, i+1) = InitVal(2,gridpt);
    end
    
    % evaluate physical flow variables from InitVal for next subinterval
    
    [Mach_D, Mrho_D, Press_D, Temp_D, Vel_D] = ...
                          Calc_FlowVarResults(Gamma,Tot_X_Pts,A,InitVal);
    
    if i == n
      % calculate relative error in primary flow variables
    
        Rel_MachErr = abs(Mach_D - Mach_E)./Mach_E;
        Rel_MrhoErr = abs(Mrho_D - Mrho_E)./Mrho_E;
        Rel_PressErr = abs(Press_D - Press_E)./Press_E;
        Rel_TempErr = abs(Temp_D - Temp_E)./Temp_E;
        Rel_VelErr = abs(Vel_D - Vel_E)./Vel_E;
        
      % calculate mean mass flow rate at final time; define array to store
      
        MeanU2 = mean(U2( :, n+1));
        
        AvU2 = zeros(1, Tot_X_Pts);
        
      % calculate mean and standard deviation of relative errors and store
      
        MeanRelTempErr = mean(Rel_TempErr);
        MeanRelMachErr = mean(Rel_MachErr);
        MeanRelMrhoErr = mean(Rel_MrhoErr);
        MeanRelPressErr = mean(Rel_PressErr);
        
        SDevRelTempErr = std(Rel_TempErr);
        SDevRelMachErr = std(Rel_MachErr);
        SDevRelMrhoErr = std(Rel_MrhoErr);
        SDevRelPressErr = std(Rel_PressErr);
        
        AvRelTempErr = zeros(1,Tot_X_Pts);
        AvRelMachErr = zeros(1,Tot_X_Pts);
        AvRelMrhoErr = zeros(1,Tot_X_Pts);
        AvRelPressErr = zeros(1,Tot_X_Pts);
        
        AvPlusSDevRelTempErr = zeros(1, Tot_X_Pts);
        AvPlusSDevRelMachErr = zeros(1, Tot_X_Pts);
        AvPlusSDevRelMrhoErr = zeros(1, Tot_X_Pts);
        AvPlusSDevRelPressErr = zeros(1, Tot_X_Pts);
        
        AvMinusSDevRelTempErr = zeros(1,Tot_X_Pts);
        AvMinusSDevRelMachErr = zeros(1,Tot_X_Pts);
        AvMinusSDevRelMrhoErr = zeros(1,Tot_X_Pts);
        AvMinusSDevRelPressErr = zeros(1,Tot_X_Pts);
        
        for col = 1:Tot_X_Pts
            AvU2(col) = MeanU2;
            
            AvRelTempErr(col) = MeanRelTempErr;
            AvRelMachErr(col) = MeanRelMachErr;
            AvRelMrhoErr(col) = MeanRelMrhoErr;
            AvRelPressErr(col) = MeanRelPressErr;
            
            AvPlusSDevRelTempErr(col) = MeanRelTempErr + ...
                                          SDevRelTempErr;
            AvPlusSDevRelMachErr(col) = MeanRelMachErr + ...
                                          SDevRelMachErr;
            AvPlusSDevRelMrhoErr(col) = MeanRelMrhoErr + ...
                                          SDevRelMrhoErr;
            AvPlusSDevRelPressErr(col) = MeanRelPressErr + ...
                                          SDevRelPressErr;
                                      
            AvMinusSDevRelTempErr(col) = MeanRelTempErr - ...
                                          SDevRelTempErr;
            AvMinusSDevRelMachErr(col) = MeanRelMachErr - ...
                                          SDevRelMachErr;
            AvMinusSDevRelMrhoErr(col) = MeanRelMrhoErr - ...
                                          SDevRelMrhoErr;
            AvMinusSDevRelPressErr(col) = MeanRelPressErr - ...
                                          SDevRelPressErr;
        end
        
    end
    
    % before returning to top of loop over subintervals i, store Taylor 
    %   polynomial for approximate solution in subinterval i in 
    %   z(d,rmaxp1, N, i)
    
    %for j = 1:N
    %    for aa = 1:d
    %        for bb = 1:rmaxp1
    %            z(aa,bb,j,i) = StoreLz(aa,bb,j);    % StoreLz picks up 
                                                    % another dimension
                                                    % associated with
                                                    % grid-points - fix
                                                    % this!
    %        end
    %    end
    %end
    
    % output initial simulation flow variables for subinterval i+1
    
    disp(['Code has completed subinterval ' int2str(i)])
    
    if i ~= n
        disp('  Initial condition for next subinterval is:')
    elseif i == n
        disp('  Final result for steady state U values are:')
    end
        
    InitVal'
    
    nextStartTime = n^(k-1)*hbar*i;
    
    disp(['Next subint start-time = ' num2str(nextStartTime)])
    
 end   

end

