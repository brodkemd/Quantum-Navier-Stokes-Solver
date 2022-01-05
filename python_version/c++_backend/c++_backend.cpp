#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <ratio>
#include <chrono>
#include <bits/stdc++.h>
#include <iostream>
#include <stdio.h>


#ifdef _WIN32
	#include <Windows.h>
#else
	#include <unistd.h>
#endif

#define TEST 0

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;

const double pi = 2*acos(0.0);
const double Tiny = pow(1*10, (-1 * 50));
const double half = ((double)1)/((double)2);

uint64_t t;

double v1, v2, magv, pofx, x, ratio, Momega;
int intx, nearx, TempTot, TotRuns;

double RANDOM(void) {
	t = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	srand(t);
	return FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN)));
}

double median(int arr[], int size){
   std::sort(arr, arr+size);
   if (size % 2 != 0) return (double)arr[size/2];
   return (double)(arr[(size-1)/2] + arr[size/2])/2.0;
}

double randQAEA(int M, double omega) {
        Momega = M * omega;

        x = -1;
        ratio = -1;

		while ((x < 0 || x > (M-1)) || (RANDOM() > ratio)) {	
			do {
				v1 = RANDOM();
				v2 = 2 * RANDOM() - 1;
				magv = pow(v1, 2) + pow(v2, 2);
			} while (magv > 1);
			
			x = v2/v1 + Momega;
			intx = round(x);
			if (intx >= 0) {
				if (intx <= (M-1)) nearx = intx;
				else nearx = M;

			}
			else { nearx = intx; }
			// pofx evaluates QAEA probability distribution at nearx
			pofx = half * pow((sin(pi*(Momega - nearx + Tiny))), (2)) / pow((M*sin((pi/M)*(Momega - nearx + Tiny))), (2)) + half*pow((sin(pi*(M - Momega - nearx + Tiny))), (2)) / pow((M*sin((pi/M)*(M - Momega - nearx + Tiny))), (2));
			
			ratio = (1 + pow((intx - Momega), (2)))*pofx;
		}

		return nearx;
}

extern "C" double QAmpEst(double M, double delta, double omega) {
	TempTot = ceil(1.25 * (-8 * log(delta)));

	if (TempTot % 2 == 0) TotRuns = TempTot + 1;
	else TotRuns = TempTot;


	int Estimates[TotRuns];
	//QAMPEST Estimates unknown quantum amplitude using QAEA.
	// start loop to carry out TotRuns simulation runs
	for (int i = 0; i < TotRuns; i++){
		Estimates[i] = randQAEA(M, omega);
	}
	
	return pow(((sin(pi*median(Estimates, TotRuns))/M)), (2));

}

int main(){
	std::cout << QAmpEst(256, 1.958004980651129e-05, 0.2) << "\n";
	return 0;
}
