#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_AR1.h"

using namespace std;

const double sigmasq = 1;
const float phi = -0.5;
const float p = 0.1;
const double N = 30;

vector < double > X;
vector < int > vect_obs_N;
vector < vector < int > > mat_B;

/* this is the code to simulate the AR(1) model and find the vector of the observations.
Sampling from a normal distribution I populate a vector X of events.
Sampling form a geometric distribution I populate a vector "vector_ti" for the times of observations
for the event that happened at time i (event x_i will be observed at time t_i from when it happened).
Then creating a matrix of observations "Obs" that will have on each row
the events that have been observed up to the time corresponding to the row number.
So, at time 0 I will have 1 element in the row that might or might not have been observed,
at time 1 I will have two elements on the row, some observed, some not, and so on
I will have 0s whenever the element have have not been yet obeserved.
The values of the parameters at this stage are fixed and are sigma^2 = 1,
phi = 0.5, p = 0.4 */


int ar1() {

	random_device rd;
	mt19937 generator(rd());

	//draw from a normal distribution to simulate the AR(1) model. Store in a vector "X"
	normal_distribution < double > normal1(0, (sigmasq / (1 - pow(phi, 2) ) ) );
	X.push_back(normal1(generator));

	for (size_t i = 1; i < N; i++) {
		normal_distribution < double > normal2( (phi * X[i - 1] ), sigmasq);
		X.push_back(normal2(generator));
	}

	
	//draw from a Bernoulli distribution to simulate our knowledge of the system
	//if 0 draw again and store results in a vector "vec_B"
	//create matrix "mat_B" with vectors "vec_B"
	//vector <  vector < int > > mat_B;
	vector < int > vec_B;
	vector < int > new_vec_B;
	bernoulli_distribution ber(p);
	int first_b = ber(generator);
	vec_B.push_back(first_b);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) {
			if (vec_B[j] == 0) {
				int b_again = ber(generator);
				new_vec_B.push_back(b_again);
			}
			else { new_vec_B.push_back(1); }
		}
		int last_b = ber(generator);
		new_vec_B.push_back(last_b);
		mat_B.push_back(new_vec_B);
		vec_B = new_vec_B;
		new_vec_B.clear();
	}
	for (const size_t i : mat_B[N - 1]) {
		vect_obs_N.push_back(i);
	}

	//for simplicity we assume we always observe the first element
	mat_B[0][0] = 1;
	


	//Creating a dat file with the values of the vector of the observed events "vect_obs_N"
	//at the current time N calling it "vector_Z.dat"
	ofstream outFile1("./vector_B.dat");
	//outFile1 << endl;
	for (int n : mat_B[N-1]) {
		outFile1 << n << endl;
	}
	outFile1.close();

	//Creating a dat file with the values of X
	//calling it "vector_X.dat""
	std::ofstream outFile2("./vector_X.dat");
	//outFile2 << endl;
	for (double n : X) {
		outFile2 << n << endl;
	}
	outFile2.close();
	
	return 0;
}
