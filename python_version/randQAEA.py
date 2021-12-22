
import numpy as np

def randQAEA(M, omega):
    #randQAEA Generates random deviate for Quantum Amplitude Estimation.
    #   Function generates random deviate for probability distribution
    #   produced by Quantum Amplitude Estimation algorithm (see Brassard
    #   et al., quant-ph/0005055.
    #
    #   returns randev in interval [0,1] M is dimension of Hilbert space
    #   omega is unknown scalar in [0,1] that determines the desired amplitude
    #   a = ( sin(pi*omega) )**2 which is calculated by QAmpEst function.
    #
    #   function based on rejection sampling described in Numerical Recipes
    #
    # ALSO: QAmpEst

    Momega = M*omega

    Tiny = 1*10**(-1 * 50)       # tiny number used to prevent 0/0 in pofx calculation

    x = -1
    ratio = -1

    # begin calculation of randev

    while ((x < 0 | x > (M-1)) | (np.random.rand(1) > ratio)):
        v1 = np.random.rand(1)
        v2 = 2*np.random.rand(1) - 1
        magv = v1**2 + v2**2

        while (magv > 1):
            v1 = np.random.rand(1)
            v2 = 2*np.random.rand(1) - 1
            magv = v1**2 + v2**2

        y = v2/v1
        x = y + Momega
        intx = np.round(x)
        inty = intx - Momega

        if (intx >= 0):
            if (intx <= (M-1)):
                nearx = intx
            else:
                nearx = M

        else:
            # 'Warning: intx is negative - set nearx to intx'

            nearx = intx


        # pofx evaluates QAEA probability distribution at nearx

        pofx = (1/2)*(np.sin(np.pi*(Momega - nearx + Tiny)))**(2)/(M*np.sin((np.pi/M)*(Momega - nearx + Tiny)))**(2)+(1/2)*(np.sin(np.pi*(M - Momega - nearx + Tiny)))**(2)/(M*np.sin((np.pi/M)*(M - Momega - nearx + Tiny)))**(2)

        ratio = (1 + (inty)**(2))*pofx

    # nearx is the desired random deviate randev

    randev = nearx

    # computation completed
    return randev
