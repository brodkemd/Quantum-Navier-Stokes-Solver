function [M, Mrho, Press, Temp, Vel, SW_JumpP, SW_JumpM] = ...
                            Calc_ExactResultspSW(Gamma, Tot_X_Pts, A, ...
                                                  Exit_Pressure, Del_x)
%CALC_EXACTRESULTSPSW assigns exact flow values when shock-wave present
%   Calc_ExactResultspSW determines the exact values of the flow variables
%       at all grid-points for 1D inviscid compressible flow through a 
%       symmetric nozzle when a shock-wave is present. All flow variables 
%       are dimensionless as in rest of simulation code. Formulas on which
%       this function is base appear in Anderson, "Modern Compressible
%       Flow", (MCF) Chapters 3 and 5. See also: Anderson, "Computational 
%       Fluid Dynamics", (CFD) Chapter 7; and entries for 07/30/19, 
%       07/31/19, and 08/01/19 in my research notebooks.
%
%   INPUTS:
%           Gamma = ratio of specific heats
%           Tot_X_Pts = number of spatial grid-points
%           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
%           Exit_Pressure = pressure at nozzle exit
%           Del_x = distance between adjacent grid-points
%
%   OUTPUTS:
%           M = 1 x Tot_X_Pts array storing Mach number at all grid-points
%           Mrho = 1 x Tot_X_Pts array storing mass density at all grd-pts
%           Press = 1 x Tot_X_Pts array storing pressure at all grd-pts
%           Temp = 1 x Tot_X_Pts array storing temperature at all grd-pts
%           Vel = 1 x Tot_X_Pts array storing velocity at all grd-pts
%           SW_JumpP = ratio of pressure behind SW/pressure ahead of SW
%           SW_JumpM = ratio of Mach number behind SW/Mach number ahead SW

% Initialize arrays to store return arrays

 M = zeros(1, Tot_X_Pts);
 Mrho = zeros(1, Tot_X_Pts);
 Press = zeros(1, Tot_X_Pts);
 Temp = zeros(1, Tot_X_Pts);
 Vel = zeros(1, Tot_X_Pts);
 
% assign values for useful parameters

 ithroat = (Tot_X_Pts + 1)/2;   % grid-point index for nozzle throat
 itp1 = ithroat + 1;   % grid-point index for start of supersonic region
 itm1 = ithroat - 1;       % grid-point index for end of subsonic region
 
 fac1 = -1/(Gamma -1);
 fac2 = (fac1)^(2);
 fac3 = 2/(Gamma - 1);
 fac4 = 2/(Gamma + 1);
 fac5 = (Gamma - 1)/2;
 
 exponent1 = (Gamma + 1)/(Gamma - 1);
 exponent2 = Gamma/(Gamma - 1);
 exponent3 = exponent1/2;
 exponent4 = 2*( 2*exponent2 - 1 );
 
 ratioPA = Exit_Pressure*A(Tot_X_Pts); % Eq.(7.126a) in CFD
 
 arg1 = fac2 + fac3*( fac4^(exponent1) )*(1/ratioPA)^(2);
 
 M(Tot_X_Pts) = sqrt( fac1 + sqrt(arg1) ); % Eq.(5.28) in MCF
 
 arg2 = 1 + fac5*( M(Tot_X_Pts) )^(2);
 
 ratioPoe = (arg2)^(exponent2);      % Eq.(7.128) in CFD
 ratioP21 = ratioPoe*Exit_Pressure;  % Eq.(7.132) in CFD
 
% throat is always at Mach 1 so ...

 M(ithroat) = 1;        % Eq.(5.20) in MCF
 
% evaluate other flow variables at throat

 argum = 1 + fac5*(M(ithroat))^(2);
 
 Press(ithroat) = (argum)^(-exponent2);     % Eq.(3.30) in MCF
 Mrho(ithroat) = (argum)^(fac1);            % Eq.(3.31) in MCF
 Temp(ithroat) = 1/argum;                   % Eq.(3.28) in MCF
 
 Vel(ithroat) = M(ithroat)*sqrt(Temp(ithroat)); % Eq. (7.76) in CFD
 
% begin calculation to locate shock-wave position
 
 term1 = (6^(exponent4))*(1/ratioP21)^(2);
 
 fcoeff7 = 855100224 - term1;
 
% assign polynomial that yields Mach number just before the shock-wave; 
%   see Eq.(7.127) in CFD. NOTE: flow is supersonic there!

 f = [ 16807 576240 8406930 67347560 316914885 fcoeff7 1123343340 ...
         214527600 -594129375 216650000 -34518750 2625000 -78125 ];
     
 z = roots(f);
 
% loop over roots: if root is real, and greater than 1 (supersonic flwo), 
%   then assign it to z1 and set M1 = sqrt(z1)

 for j = 1: exponent4
     if isreal(z(j))
         if z(j) > 1
             z1 = z(j);
             M1 = sqrt(z1);
         end
     end
 end
 
 % calculate jump in pressure and Mach number across shock-wave
 
 SW_JumpP = 1 + ((2*Gamma)/(Gamma+1))*( M1^2 - 1 );  % Eq.(7.137) in CFD
 SW_JumpM = sqrt( (1 + fac5*M1^2)/(Gamma*M1^2 - fac5) ); % Eq.(7.138) CFD

 arg3 = 1 + fac5*z1;
 arg4 = fac4*arg3;
 
% calculate nozzle area Asw at shock-wave location

 Asw = (1/M1)*( arg4 )^(exponent3);     % Eq.(5.20) in MCF
 
% calculate position Xsw of shock-wave - use root greater than 1.5 since
%   shock-wave is in divergent part of nozzle

 Xsw = 1.5 + sqrt( (Asw - 1)/2.2 );     % solve Eq.(7.73) in CFD for x
 
% calculate grid-point index isw of shock-wave location

 isw = round( 1 + Xsw/Del_x);
 
% define grid-point label for point just downstream of shock-wave
 
 iswp1 = isw + 1;
 
% now that we know where the shock-wave is we can determine the exact
%   values of the flow variables at all grid-points.

% begin by calculating flow variables upstream of shock-wave. First
%   need to determine Mach number at all upstream grid-points

% flow is subsonic before throat

 for i = 1:itm1
     term2 = (6^(exponent1))*(A(i)^(2));
     
     gcoeff1 = 18750 - term2;
     
     g = [ 1 30 375 2500 9375 gcoeff1 15625 ];
     
     z = roots(g);
     
     % want root that's real, positive, and less than 1 (subsonic flow)
     
     for j = 1:exponent1
         if isreal(z(j))
             if z(j) > 0
                 if z(j) < 1
                     M(i) = sqrt(z(j));
                 end
             end
         end
     end
 end
     
% now repeat for supersonic flow upstream of shock-wave

 for i = itp1: isw
     term2 = ( 6^(exponent1) )*( A(i)^(2) );
     
     gcoeff1 = 18750 - term2;
     
     g = [ 1 30 375 2500 9375 gcoeff1 15625 ];
     
     z = roots(g);
     
     % want root that's real and greater than 1 (supersonic flow)
     
     for j = 1:exponent1
         if isreal(z(j))
             if z(j) > 1
                 M(i) = sqrt(z(j));
             end
         end
     end
 end
 
% now calculate flow variables upstream of shock-wave

 for i = 1:isw
     arg5 = 1 + fac5*(M(i))^(2);
     
     Press(i) = (arg5)^(-exponent2);
     Mrho(i) = (arg5)^(fac1);
     Temp(i) = 1/arg5;
     
     Vel(i) = M(i)*sqrt(Temp(i));
 end
 
% now calculate flow downstream of shock-wave; flow is subsonic here

% start by calculating sonic area AstarE for nozzle exit

 arg7 = fac4*( 1 + fac5*( M(Tot_X_Pts) )^(2) );
 
% derivation of following formula in my research notebook, 07-31-19 entry,
%   page 11
 
 AstarE = A(Tot_X_Pts)*M(Tot_X_Pts)*( arg7 )^(-exponent3);
 
% begin loop over downstream grid-points

for i = iswp1:Tot_X_Pts
    term3 = ( (6^(exponent1))*(A(i))^(2) )/( AstarE^(2) );
    
    hcoeff1 = 18750 - term3;
    
    h = [ 1 30 375 2500 9375 hcoeff1 15625 ];
    
    z = roots(h);
    
    % want root that's real, positive , and less than 1 (subsonic flow)
    
    for j = 1:exponent1
        if isreal(z(j))
            if z(j) > 0
                if z(j) < 1
                    M(i) = sqrt(z(j));
                end
            end
        end
    end
end  
    
% calculate flow variables downstream of shock-wave

 for i = iswp1:Tot_X_Pts
     arg6 = 1 + fac5*(M(i))^(2);
     
     % following formulas derived in my research notebook, 07-30-19 entry,
     %      page 7, Eqs.(15a-c)
     
     Press(i) = ratioP21*(arg6)^(-exponent2);
     Mrho(i) = ratioP21*(arg6)^(fac1);
     Temp(i) = 1/arg6;
     
     Vel(i) = M(i)*sqrt(Temp(i));
 end
 
 % finally, rescale Mrho, Press, and Temp by their values at nozzle 
 %      entrance so boundary values there are unity as in simulation BCs
 %{
 Mrho = Mrho/Mrho(1);
 Press = Press/Press(1);
 Temp = Temp/Temp(1);
 Vel = M.*sqrt(Temp);
 %}
end

