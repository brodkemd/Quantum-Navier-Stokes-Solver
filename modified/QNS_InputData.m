% An explanation of the physical meaning of the variables below is given 
%  in the comment section at the beginning of the function 
%  QNavStokes_Solvr.m. 

a = 0;
Tot_TSteps = 1400;
in_n = 16;
d = 3;
r = 2;
err1 = 0.005;
delta = 0.005;
rho = 1;
x_min = 0;
x_max =3.0;
Tot_X_Pts = 61;
Gamma = 1.4;
Exit_Pressure = 0.6784;
Shock_Flag = 1;

% In_Mass_Flow set to exact value for sub-to-supersonic flow without
%   shock-wave

In_Mass_Flow = 0.579;

% ICMFlowErrScale less than or about 0.01 which is size of shift used by
%   Anderson

ICMFlowErrScale = 0.00;

% ICrhoErrScale is 10% (2%) of minimum density in steady-state solution
%   when shockwave absent (present)

ICrhoErrScale = 0.02;

% ICtempErrScale is 2% (1%) of minimum temperature in steady-state 
%   solution when shockwave absent (present)

ICtempErrScale = 0.01;


% settings for the QAmpEst algorithm 
server_compute = false; % computes the value on the server
server_log = true;

for i = 1:16
    in_n = i;
    QNavStokes_Solvr( a, Tot_TSteps, in_n, d, r,err1, ...
                                        delta, rho, x_min, x_max, ...
                                        Tot_X_Pts, Gamma, ...
                                        ICMFlowErrScale, ICrhoErrScale, ...
                                        ICtempErrScale, ...
                                        In_Mass_Flow, ...
                                        Exit_Pressure, Shock_Flag, server_compute, server_log)
end

