#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_AR1.h"
#include "numeric"

using namespace std;

void F_print_matrix(vector < vector < int > > m);
void F_print_matrix(vector < vector < double > > M);
void F_print_vector(vector < double > v);
void F_print_vector(vector < int > V);

//This is the code for the method:
//Firstly I calculate the lower triangular 3-dimentional matrices called for the sampled events
//and for the samples with substitutions for every particle.
//I then calculate the lower triangular 3-dimensional matrices for the sampled observations
//and for the sampled observations with substitutions
//Finally I calculate the weights.


int main() {

	random_device rd;
	mt19937 generator(rd());

	ar1();


	//DEFINITIONS

	//define the number of particles
	int n = 1000;
	//define the container for the sampled events and the sampled observations (0s and 1s)
	vector < vector < vector < double > > > sample(N, vector < vector < double > > (N, vector < double > (n, 0.0)));
	//define the container for the new sampled events and the new sampled observations (0s and 1s)
	vector < vector < vector < double > > > corr_sample(N, vector < vector < double > >(N, vector < double > (n, 0.0)));
	vector < vector < vector < double > > > resampled(N, vector < vector < double > >(N, vector < double > (n, 0.0)));
	//define and initialise the containter for the unnormalised weights
	vector < vector < double > > un_weights(n, vector < double > (N, 0.0));
	//define and initialise the container for the normalised weights
	vector < vector < double > > weights(n, vector < double >(N, 0.0));

	
	
	//initialize:
	//draw from a normal distribution with mean 0
	//and variance sigma^2/(1-phi^2) for every particle
	//here we make the assumption that the first nest is never observed
	//therefore I don't need to make a correction at initialization
	for (int i = 0; i < n; i++) {
		normal_distribution < double > normalDist(0, sigmasq / (1 - phi * phi));
		un_weights[i][0] = 1;
		weights[i][0] = 1.0 / n;
		sample[0][0][i] = normalDist(generator);
		corr_sample[0][0][i] = sample[0][0][i];
	}
 
	//Iterate	
	for (int j = 1; j < N; j++) {
		for (int i = 0; i < n; i++) {
			vector < double > sum_vec;
			//draw the next value for x from the transition distribution
			normal_distribution < double > normalDist(phi * (corr_sample[j][j-1][i]), sigmasq);
			sample[j][j][i] = normalDist(generator);
			corr_sample[j][j][i] = sample[j][j][i];
			for (int k = 1; k < j + 1; k++) {
				sample[j][k - 1][i] = corr_sample[j - 1][k - 1][i];
				//make the corrections
				if (mat_B[j][k - 1] == 1) {corr_sample[j][k - 1][i] = X[k - 1]; }
				else { corr_sample[j][k - 1][i] = corr_sample[j - 1][k - 1][i]; }
				//calculate the weights
				if ( ( (mat_B[j - 1][k - 1] == mat_B[j][k - 1]) ) ) {
					double num_arg = pow( (corr_sample[j][k][i] - phi * corr_sample[j][k-1][i]), 2);
					double den_arg = pow( (sample[j][k][i] - phi * sample[j][k - 1][i]), 2);
					double log_elem = den_arg - num_arg;
					sum_vec.push_back(log_elem);
				}
			}
			double sum_of_logs = accumulate(sum_vec.begin(), sum_vec.end(), 0.0);
			double W = exp(sum_of_logs);
			un_weights[i][j] = W;
		} 
		//normalise the weights
		double sum_of_weights{ 0 };
		for (int l = 0; l < n; l++) {
			sum_of_weights = sum_of_weights + un_weights[l][j];
		}
		for (int i = 0; i < n; i++) {
			weights[i][j] = un_weights[i][j] / sum_of_weights;
		}
		//resampling
		vector < double > drawing_vector(n, 0.0);
		for (int i = 0; i < n; i++) {
			drawing_vector[i] = weights[i][j];
		}
		//F_print_vector(drawing_vector);
		for (int k = 0; k < j + 1; k++){
			for (int i = 0; i < n; i++) {
				discrete_distribution < int > discrete(drawing_vector.begin(), drawing_vector.end());
				resampled[j][k][i] = corr_sample[j][k][discrete(generator)];
			}
		}
	} 

	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			cout << sample[29][i][j] << " ";
		}
		cout << endl;
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			cout << corr_sample[29][i][j] << " ";
		}
		cout << endl;
	}
	cout << "matrix B" << endl;
	F_print_matrix(mat_B);


	//Create a .csv file with the resampled particles (transposing)
	ofstream outFile("./resampled_000.csv");
		outFile << endl;
		vector < vector < double > > tresampled(N, vector < double >(n, 0.0));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < N; j++) {
				tresampled[j][i] = resampled[N - 1][j][i];
				outFile << tresampled[j][i] << ",";
			}
		outFile << endl;
	}
	outFile.close();

	//create a .csv file containing the parameters
	ofstream outparam("./parameters.csv");
	outparam << "sigmasq" << "," << "phi" << "," << "p" << "," << "N" << "," << "part" << endl;
	outparam << sigmasq << "," << phi << "," << p << "," << N << "," << n << "," << endl;
	outparam.close();


	gold();

	return 0;



}






