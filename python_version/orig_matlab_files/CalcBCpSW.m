function U_Bvals = CalcBCpSW(U,A,Gamma,d,Tot_X_Pts,Exit_Pressure)
%CALCBCPSW evaluates flow variables at nozzle boundaries-shock-wave present
%   CalcBCpSW evaluates flow variables at boundaries of a nozzle. Flow is
%       1D compressible, inviscid flow with shock-wave present.
%   
%   INPUTS:
%           U = d x Tot_X_Pts array with current values of flow variables
%                   at interior grid-points only. Values at boundary
%                   grid-points are zero - will be updated in calling
%                   program.
%           A = 1 x Tot_X_Pts array with values of nozzle area at
%                   grid-points.
%           Gamma = ratio of specific heats
%           d = number of flow variables at a grid-point
%           Tot_X_Pts = number of grid-points
%           Exit_Pressure = pressure at nozzle exit
%
%   OUTPUT:
%           U_Bvals = d x 2 array with values of flow variable at nozzle
%                       entrance and exit.

U_Bvals = zeros(d,2);

fac1 = 1/(Gamma -1);
fac2 = Gamma/2;

GP_Nm1 = Tot_X_Pts -1;
GP_Nm2 = Tot_X_Pts -2;

% evaluate U at nozzle entrance:

U_Bvals(1,1) = A(1);
U_Bvals(2,1) = 2*U(2,2) - U(2,3);

fac3 = (U_Bvals(2,1)/U_Bvals(1,1))^(2);

U_Bvals(3,1) = U_Bvals(1,1)*(fac1 + fac2*fac3);

% evaluate flow at nozzle exit:

for k = 1:(d-1)
    U_Bvals(k,2) = 2*U(k,GP_Nm1) - U(k,GP_Nm2);
end

U_Bvals(3,2) = Exit_Pressure*A(Tot_X_Pts)*fac1 ...
                    +( fac2*( U_Bvals(2,2))^(2) )/U_Bvals(1,2);

end

