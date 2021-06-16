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

	//this function is the realisation of the AR(1) model
	//it is generated in AR1_model.cpp
	ar1();


	//DEFINITIONS

	//define the number of particles
	int n = 3000;
	//define the container for the sampled events and the sampled observations (0s and 1s)
	vector < vector < vector < double > > > sample(N, vector < vector < double > > (N, vector < double > (n, 0.0)));
	//define the container for the new sampled events and the new sampled observations (0s and 1s)
	vector < vector < vector < double > > > corr_sample(N, vector < vector < double > > (N, vector < double > (n, 0.0)));
	vector < vector < vector < double > > > resampled(N, vector < vector < double > > (N, vector < double > (n, 0.0)));
	//define and initialise the containter for the unnormalised weights
	vector < vector < double > > un_weights(n, vector < double > (N, 0.0));
	//define and initialise the container for the normalised weights
	vector < vector < double > > weights(n, vector < double > (N, 0.0));

	
	
	//initialize:
	//draw from a normal distribution with mean 0
	//and variance sigma^2/(1-phi^2) for every particle
	//here we make the assumption that the first nest is never observed
	//therefore I don't need to make a correction at initialization
	for (int i = 0; i < n; i++) {
		un_weights[i][0] = 1;
		weights[i][0] = 1.0 / n;
		for (int j = 0; j < N; j++){	
			sample[j][0][i] = X[0];
			resampled[j][0][i] = X[0];
			corr_sample[j][0][i] = X[0];
		}
	}
 
	//Iterate	
	for (int j = 1; j < N; j++) {
		for (int i = 0; i < n; i++) {
			vector < double > sum_vec;
			for (int k = 0; k < j; k++) { 
				sample[j][k][i] = resampled[j - 1][k][i];
				corr_sample[j][k][i] = resampled[j - 1][k][i];
			}
			normal_distribution < double > normalDist1((phi * (sample[j][j-1][i])), sigmasq);
			sample[j][j][i] = normalDist1(generator);
			corr_sample[j][j][i] = sample[j][j][i];
			//make the corrections
			for (int k = 0; k < j + 1; k++) {
				if (mat_Z[j][k] == 1) { corr_sample[j][k][i] = X[k]; }
			}
			//calculate the weights
			//this condition ensures that we are only calculating the partial weights that are not 1
			for (int k = 1; k < j + 1; k++) {
				if ( (corr_sample[j][k][i] != sample[j][k][i]) || (corr_sample[j][k - 1][i] != sample[j][k - 1][i]) ) {
					double num_arg = pow((corr_sample[j][k][i] - phi * corr_sample[j][k - 1][i]), 2);
					double den_arg = pow((sample[j][k][i] - phi * sample[j][k - 1][i]), 2);
					double log_elem = (1 / (2 * sigmasq)) * (num_arg - den_arg);
					sum_vec.push_back(log_elem);
				}
				else { sum_vec.push_back(0); }
			}			
			double sum_of_logs = accumulate(sum_vec.begin(), sum_vec.end(), 0.0);
			double W = exp(sum_of_logs);
			un_weights[i][j] = W;
		} 
		//normalise the weights
		double sum_of_weights{ 0 };
		for (int i = 0; i < n; i++) {
			sum_of_weights = sum_of_weights + un_weights[i][j];
		}
		for (int i = 0; i < n; i++) {
			weights[i][j] = un_weights[i][j] / sum_of_weights;
		}
		//resampling
		vector < double > drawing_vector(n, 0.0);
		for (int i = 0; i < n; i++) {
			drawing_vector[i] = weights[i][j];
		}
		F_print_vector(drawing_vector);
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < j + 1; k++) {
				resampled[j][k][i] = corr_sample[j][k][i];
			}
		}
		for (int i = 0; i < n; i++) {
			discrete_distribution < int > discrete(drawing_vector.begin(), drawing_vector.end());
			resampled[j][j][discrete(generator)];
		}
	}
	/*
	cout << "sample" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			cout << sample[29][i][j] << " ";
		}
		cout << endl;
	}
	cout << "corr_sample 1" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			cout << corr_sample[29][i][j] << " ";
		}
		cout << endl;
	}
	cout << "corr_sample 2" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			cout << corr_sample[28][i][j] << " ";
		}
		cout << endl;
	}
	cout << "resampled" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			cout << resampled[29][i][j] << " ";
		}
		cout << endl;
	}
	cout << "matrix B" << endl;
	F_print_matrix(mat_Z);
	*/

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

	cout << "end of corrections" << endl;

	//this function is the gold standard for the AR(1) model
	//it is generated in AR1_sim_and_int.cpp
	gold();

	cout << "end of gold" << endl;

	return 0;



}






