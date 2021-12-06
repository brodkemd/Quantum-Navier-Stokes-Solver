function [tt, hhbar] = IPrtn( Init, Fnl, nn, NN )
%IPRTN partitions interval [Init, Fnl] into nn*NN subintervals
%   IPrtn partitions [Init, Fnl] into nn*NN subintervals by introducing
%   NN - 1 intermediate times tt(j,i). It does a primary partition
%   of [Init, Fnl] into nn subintervals and then partitions each primary
%   subinterval into NN sub-subintervals. Primary subintervals are
%   indexed by 1 <= i <= nn, and the sub-subintervals are indexed by
%   1 <= j <= NN. The time tt(j,i) is the rightmost time in
%   sub-subinterval j in primary subinterval i. (NOTE: NN = nn^(k-1),
%   where k is the number of recursion levels for algorithm A_k.)

hh = (Fnl - Init)/nn;     % width of a primary subinterval
hhbar = hh/NN;             % width of a sub-subinterval

tt = zeros(NN,nn);        % create/initialize times array; first (second)
                        % slot specifies sub- (primary) subinterval

% fill in tt array-element values; Note tt(N,n) = Fnl

for i = 1:nn
    for j = 1:NN
        if i == 1
            tt(j,i) = Init + j*hhbar;    % fill times in first subinterval
        else
            tt(j,i) = tt(NN,i-1) + j*hhbar; % fill remaining subintervals
        end
    end
end

end

