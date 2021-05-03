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
const float phi = 0.5;
const float p = 0.1;
const double N = 30;

vector < double > X;
vector < vector < size_t > > obs;
vector < size_t > vect_obs_N;

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

	//Sampling from a normal distribution simulate the AR(1) model. Put the values in a vector "X"
	normal_distribution < double > normalDist(0, sigmasq / (1 - phi * phi));
	X.push_back(normalDist(generator));

	for (size_t i = 1; i < N; i++) {
		normal_distribution < double > normalDist(phi * X[i - 1], sigmasq);
		X.push_back(normalDist(generator));
	}

	//Sampling from a geometric distribution to get the values of the number of Bernoulli trials
	//needed to get an observation for each event.
	//Put the values in a vector called "vector_ti".
	vector < double > vector_ti;
	for (size_t i = 0; i < N; i++) {
		geometric_distribution <> geoDist(p);
		size_t ti = geoDist(generator);
		vector_ti.push_back(ti);
	}

	//Populating the matrix of observations of "0" and "1" callled "obs"
	//using the previusly calculated "vector_ti".
	//Each line of the vector is the state of the observations at the corresponing time.
	//Populate the vector "vect_obs_N" for the final time
	for (size_t j = 0; j < N; j++) {
		vector < size_t > tempvec;
		for (size_t i = 0; i < j + 1; i++) {
			if (vector_ti[i] <= j - i) {
				tempvec.push_back(1);
			}
			else {
				tempvec.push_back(0);
			}
		}
		obs.push_back(tempvec);
	}
	vect_obs_N = obs[N - 1];

	//Creating a dat file with the values of the vector of the observed events "vect_obs_N"
	//at the current time N calling it "real_data.dat"
	ofstream outFile1("./real_data.dat");
	outFile1 << endl;
	for (double n : vect_obs_N) {
		outFile1 << n << endl;
	}
	outFile1.close();

	//Creating a dat file with the values of X
	//calling it "vector_X.dat""
	std::ofstream outFile2("./vector_X.dat");
	outFile2 << endl;
	for (double n : X) {
		outFile2 << n << endl;
	}
	outFile2.close();

	return 0;
}