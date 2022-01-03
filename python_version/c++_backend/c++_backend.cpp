#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <ratio>
#include <chrono>
#include <bits/stdc++.h>
#include<iostream>

#define TEST 1

using namespace std;

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;

const double pi = 2*acos(0.0);
double Tiny = 1*10^(-1 * 50);

double v1, v2, magv, y, nearx, pofx, ratio;
int intx, inty;

uint64_t t;

double RANDOM(void) {
	t = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

	//std::cout << t << "\n";
	srand((unsigned) t);
	return FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN)));
}

double median(int arr[], int size){
   sort(arr, arr+size);
   if (size % 2 != 0)
      return (double)arr[size/2];
   return (double)(arr[(size-1)/2] + arr[size/2])/2.0;
}

double randQAEA(double M, double omega) {
		//srand((unsigned) time(0));
        double Momega = M*omega;

        double x = -1;
        double ratio = -1;
		
		while ((x < 0 || x > (M-1)) || (RANDOM() > ratio)) {
			cout << "looping\n";
			v1 = RANDOM();
			v2 = 2*RANDOM() - 1;
			magv = pow(v1, 2) + pow(v2, 2);
			
			while (magv > 1) {
				v1 = RANDOM();
				v2 = 2*RANDOM() - 1;
				magv = pow(v1, 2) + pow(v2, 2);
			}
			
			y = v2/v1;
			x = y + Momega;
			intx = round(x);
			inty = intx - Momega;
			
			if (intx >= 0) {
				if (intx <= (M-1)) {
					nearx = intx;
				}
				else {
					nearx = M;
				}
			}
			else {
			// 'Warning: intx is negative - set nearx to intx'
				nearx = intx;
			}
			
			// pofx evaluates QAEA probability distribution at nearx
			
			pofx = (1/2)*pow((sin(pi*(Momega - nearx + Tiny))), (2))/
            pow((M*sin((pi/M)*(Momega - nearx + Tiny))), (2))
            + (1/2)*pow((sin(pi*(M - Momega - nearx + Tiny))), (2)) /
            pow((M*sin((pi/M)*(M - Momega - nearx + Tiny))), (2));

			ratio = (1 + pow((inty), (2)))*pofx;
		}

		return nearx;
}
extern "C" double QAmpEst(double M, double delta, double omega) {
#if TEST
	int TempTot = ceil(1.25 * (-8 * log(delta)));
	int TotRuns;

	if (TempTot % 2 == 0) TotRuns = TempTot + 1;
	else TotRuns = TempTot;

	cout << TotRuns << "\n";
	
	int Estimates[TotRuns];
	//QAMPEST Estimates unknown quantum amplitude using QAEA.
	// start loop to carry out TotRuns simulation runs
	for (int i = 0; i < TotRuns; i++){
		Estimates[i] = randQAEA(M, omega);
		cout << i << endl;
	}
	
	return pow(((sin(pi*median(Estimates, TotRuns))/M)), (2));
#else
	return RANDOM();
#endif

}

int main(){
	std::cout << QAmpEst(256, 1.958004980651129e-05, 0.2) << "\n";

	return 0;
}
