function [aEstimate, trueValue, error, message] =...
                              QAmpEst(M, delta, omega, server)
%QAMPEST Estimates unknown quantum amplitude using QAEA.
%   Function uses Quantum Amplitude Estimation algorithm (QAEA)
%   to estimate the unknown quantum amplitude a = (sin(pi*omega))^2.
%   (See Brassard et al., quant-ph/0005055).
%
%   (i) Function estimates number of runs (TotRuns) needed to insure
%       error probability for final estimate is less than delta.
%       TotRuns calculated in two steps: (a) TempTot = 
%       ceil(FudgeFactor*(-8*log(delta))), where FudgeFactor reduces
%       probability that an estimate won't satisfy the QAEA upper bound;
%       (b) TotRuns equal smallest odd integer greater than TempTot which
%       insures median calculation returns an element of the Estimates
%       array.
%   (ii) Function generates TotRun estimates for y (=M*omega/pi) and
%        stores them in Estimates[TotRuns] array. Median of this array
%        is estimate of y; amplitude estimate aEstimate = (sin(pi*y/M))^2.
%   (iii) Error on final estimate is O(1/M).
%
%   ALSO: randQAEA.

  % sending to the server
  write(server, single(omega));

  % getting the data back from the server
  %A = read(server, 1, 'single');

  % the estimated value
  %aEstimate = (sin(pi*A))^(2);

% true value of unknown amplitude

  a = (sin(pi*omega))^2;
  
  trueValue = a;
  
% calculate total number of amplitude estimates needed (TotRuns)

  FudgeFactor = 1.25;

  TempTot = ceil(FudgeFactor*(-8*log(delta)));
  
  if (mod(TempTot,2) == 0)
      TotRuns = TempTot + 1;
  else
      TotRuns = TempTot;
  end
  
  Estimates = zeros(1,TotRuns);
  
  % start loop to carry out TotRuns simulation runs
  
  for runs = 1:TotRuns
      randev = randQAEA(M,omega);    % randQAEA generates random deviate 
                                     % with probability distribution
                                     % produced by QAEA
      
      Estimates(runs) = randev;
  end
  
  EstimateMedian = median(Estimates);

  aEstimate = (sin(pi*EstimateMedian/M))^(2);

  error = abs(aEstimate - trueValue);
  
end

