import math
import numpy as np
import randQAEA

def QAmpEst(M, delta, omega):
    #QAMPEST Estimates unknown quantum amplitude using QAEA.
    #   Function uses Quantum Amplitude Estimation algorithm (QAEA)
    #   to estimate the unknown quantum amplitude a = (sin(pi*omega))^2.
    #   (See Brassard et al., quant-ph/0005055).
    #
    #   (i) Function estimates number of runs (TotRuns) needed to insure
    #       error probability for final estimate is less than delta.
    #       TotRuns calculated in two steps: (a) TempTot = 
    #       ceil(FudgeFactor*(-8*log(delta))), where FudgeFactor reduces
    #       probability that an estimate won't satisfy the QAEA upper bound
    #       (b) TotRuns equal smallest odd integer greater than TempTot which
    #       insures median calculation returns an element of the Estimates
    #       array.
    #   (ii) Function generates TotRun estimates for y (=M*omega/pi) and
    #        stores them in Estimates[TotRuns] array. Median of this array
    #        is estimate of y amplitude estimate aEstimate = (sin(pi*y/M))^2.
    #   (iii) Error on final estimate is O(1/M).
    #
    #   ALSO: randQAEA.

    # true value of unknown amplitude
    
    a = (np.sin(np.pi*omega))**2

    trueValue = a

    # calculate upper bound on estimate error

    UpprBnd = (2*np.pi/M)*np.sqrt(a*(1-a)) + (np.pi/M)**2

    # calculate total number of amplitude estimates needed (TotRuns)

    FudgeFactor = 1.25

    TempTot = np.ceil(FudgeFactor*(-8*np.log(delta)))

    if (TempTot % 2 == 0):
        TotRuns = TempTot + 1
    else:
        TotRuns = TempTot

    Estimates = np.zeros((TotRuns))

    # start loop to carry out TotRuns simulation runs

    for runs in range(TotRuns):
        randev = randQAEA(M,omega)    # randQAEA generates random deviate 
                                    # with probability distribution
                                    # produced by QAEA

        Estimates[runs] = randev

    EstimateMedian = np.median(Estimates)

    aEstimate = (np.sin(np.pi*EstimateMedian/M))**(2)

    error = np.abs(aEstimate - trueValue)

    # no longer interested in this check
    
    if error <= UpprBnd:
        message = 'Estimate error upper bound satisfied!'
    elif error > UpprBnd:
        message = 'Estimate error upper bound violated! '
    
    return aEstimate, trueValue, error, message