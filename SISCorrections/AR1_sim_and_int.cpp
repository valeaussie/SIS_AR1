#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_AR1.h"

using namespace std;

void F_print_matrix(vector < vector < size_t > > m);
void F_print_matrix(vector < vector < double > > M);
void F_print_vector(vector < double > v);
void F_print_vector(vector < size_t > v);
void F_outExp(vector <double> exp);
void F_outVar(vector <double> var);
double F_den(int m, int tau, double phi);



int gold() {

	random_device rd;
	mt19937 generator(rd());

	/*
	//new code to run gold standard only wiht cin
	const double sigmasq = 1;
	const float phi = 0.55;
	const float p = 0.1;
	int N = 30;
	vector <size_t> vect_obs_N;
	vector <double> X;

	//opening file with simualted observations
	ifstream fin1("real_data_01_025.dat");
	if (!fin1) {
		cout << "unable to open file";
		exit(1); //terminate with error
	}
	size_t num1;
	while (fin1 >> num1)
		vect_obs_N.push_back(num1);

	//opening file with simualted X
	ifstream fin2("vector_X_01_025.dat");
	if (!fin2) {
		cout << "unable to open file";
		exit(1); //terminate with error
	}
	double num2;
	while (fin2 >> num2) X.push_back(num2);


	//end of new code to run gold standard only wiht cin
	*/

	/* This is the code to find the analytic estimates of expectation and variance
	   of the missing values in the AR(1) model. These estimates are used as
	   gold standard to evaluate our method.*/
	

	//vectors definition
	vector < size_t > missing{};
	vector < double > expectations{ X[0] };
	vector < double > variances{ 0 };


	//case 1: One observation at the beginning One at the end
	if ( (vect_obs_N[0] == 1) && (vect_obs_N[N-1] == 1) ) {
		cout << "case 1" << endl;
		for ( size_t i = 1; i < N; i++ ){
			if ( (vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 1)) {
			expectations.push_back(X[i]);
			variances.push_back(0);
			}
			else if ((vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 0)){}
			else {
				//the missing vector will contain two elements:
				//the last observed, and the next observed.
				missing.push_back(i);
				if (missing.size() == 2) {
				size_t tau = missing[0];
				size_t m = missing[1];

				//calculating the denomiator for variances and expectations
				double den = F_den(m, tau, phi);

				//calculating the numerator of expectations and variances
				for (size_t j = tau; j < m; j++) {
					//calculating the first summation for the numerator of the expectation
					//this is also a term in the numeratior of the variance
					vector < double > vecsum1{};
					for (size_t k = 0; k < (j - tau + 1); k++) {
						double power = pow(phi, 2 * k);
						vecsum1.push_back(power);
					}
					double sum1{};
					for (auto& n : vecsum1) sum1 += n;
					//calculating the power of phi for the first term in the nominator of the expectation
					double power1 = pow(phi, m - j);
					//calculating the first term in the nominator of the expectation
					double numexp1 = power1 * sum1 * X[m];
					vecsum1.clear();

					//calculating the second summation for the numerator of the expectation
					//this is also a term in the numeratior of the variance
					vector < double > vecsum2{};
					for (size_t k = 0; k < (m - j); k++) {
						double power = pow(phi, 2 * k);
						vecsum2.push_back(power);
					}
					double sum2{};
					for (auto& n : vecsum2) sum2 += n;
					//calculating the power of phi for the second term in the nominator of the expectation
					double power2 = pow(phi, j - tau);
					//calculating the first term in the nominator of the expectation
					double numexp2 = power2 * sum2 * X[tau];
					vecsum2.clear();

					//calculating the expectation
					double expect = (numexp1 + numexp2) / den;
					expectations.push_back(expect);
					//calculating the variance
					double variance = sigmasq * sum1*sum2 / den;
					variances.push_back(variance);
				}
				expectations.push_back(X[m]);
				variances.push_back(0);
				missing.clear();
				}
				else { continue; }
			
		    }
		}
		//Create a dat files for expectations and variances
		F_outExp(expectations);
		F_outVar(variances);
	}


	
	//case 2: One observation at the beginning None at the end
	else if ( (vect_obs_N[0] == 1) && (vect_obs_N[N-1] == 0) ) {
		cout << "case 2" << endl;
		size_t lastvalue{};
		for (size_t i = N-1; i >= 0; --i){
			if ( vect_obs_N[i] == 1 ) {
				lastvalue = i+1;
				break;
		    }
		}
		for ( size_t i = 1; i < lastvalue; i++ ){
			if ( (vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 1)) {
				expectations.push_back(X[i]);
				variances.push_back(0);
			}
			else if ((vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 0)){}
			else {
				//the missing vector will contain two elements:
				//the last observed, and the next observed.
				missing.push_back(i);

				//calculations for variances and expectations of middle points
				if (missing.size() == 2) {
					size_t tau = missing[0];
					size_t m = missing[1];

					//calculating the denomiator for variances and expectationsor
					double den = F_den(m, tau, phi);

					//calculating the numerator of expectations and variances
					for (size_t j = tau; j < m; j++){
						//calculating the first summation for the numerator of the expectation
						//this is also a term in the numeratior of the variance
						vector < double > vecsum1{};
						for (size_t k = 0; k < (j - tau + 1); k++) {
							double power = pow (phi, 2*k);
							vecsum1.push_back(power);
						}
						double sum1{};
						for (auto& n : vecsum1) sum1 += n;		
						//calculating the power of phi for the first term in the nominator of the expectation
						double power1 = pow(phi, m - j);
						//calculating the first term in the nominator of the expectation
						double numexp1 = power1 * sum1 * X[m];
						vecsum1.clear();

						//calculating the second summation for the numerator of the expectation
						//this is also a term in the numeratior of the variance
						vector < double > vecsum2{};
						for (size_t k = 0; k < (m - j); k++) {
							double power = pow(phi, 2*k);
							vecsum2.push_back(power);
						}
						double sum2{};
						for (auto& n : vecsum2) sum2 += n;
						//calculating the power of phi for the second term in the nominator of the expectation
						double power2 = pow(phi, j - tau);
						//calculating the first term in the nominator of the expectation
						double numexp2 = power2 * sum2 * X[tau];
						vecsum2.clear();

						//calculating the expectation
						double expect = (numexp1 + numexp2) / den;
						expectations.push_back(expect);
						//calculating the variance
						double variance = sigmasq*sum1*sum2 / den;
						variances.push_back(variance);
					}			
					expectations.push_back(X[m]);
					variances.push_back(0);
					missing.clear();
				}
			else {continue;}
			}
		}
		//estimation of the last missing points
		vector < double > vecsum{};
		for (size_t i = 0; i < N-lastvalue; i++) {
			double expect{};
			double variance{};
			double power{};
			expect = pow(phi,(i+1)) * X[lastvalue-1];
			power = pow(phi, 2 * i);
			expectations.push_back(expect);
			vecsum.push_back(power);
			double sum{};
			for (auto& n : vecsum) sum += n;
			variance = sigmasq*sum;
			variances.push_back(variance);
		}
		//Create a dat files for expectations and variances
		F_outExp(expectations);
		F_outVar(variances);
	}

	   
	/****************case 3: No observation at the beginning One at the end ************/
	else if ( (vect_obs_N[0] == 0) && (vect_obs_N[N-1] == 1) ) {
		 // these are estimations of the first missing point. Step needed because
		 // I don't know the value of the first time
		 size_t firstvalue{};
		 for ( size_t i = 0; i < N; i++ ){
		   if ( vect_obs_N[i] == 1 ) {
		 firstvalue = i;
		 break;
		   }
		 }
		 for ( size_t i = 0; i < firstvalue; i++ ){
		   double expect{};
		   double variance{};
		   double power{};
		   vector < double > vecsum{};
		   expect = pow(phi,(firstvalue-i)) * X[firstvalue];
		   power = pow(phi,2*(firstvalue-i-1));
		   expectations.push_back(expect);
		   vecsum.push_back(power);
		   double sum{};
		   for (auto& n : vecsum)
		 sum += n;
		   variance = sum*sigmasq;
		   variances.push_back(variance);
		 }
		 // end of the estimations of the first missing points.
		 expectations.push_back(X[firstvalue]);
		 for ( size_t i = firstvalue+1; i < N; i++ ){
		   if ( (vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 1)) {
		 expectations.push_back(X[i]);
		 variances.push_back(0);
		   }
		   else if ((vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 0)){}
		   else {
			 missing.push_back(i);
		 //these are all the calculations for variances and expectations
		 if (missing.size() == 2) {
		   size_t tau = missing[0];
		   size_t m = missing[1];
		   //beginning calculating the denomiator for variances and expectations
		   double den = F_den(m, tau, phi);
		   //ending calculating the denomiator for variances and expectations
		   double power1{};
		   double power2{};
		   vector < double > vecsum1{};
		   vector < double > vecsum2{};
		   double numexp1{};
		   double numexp2{};
		   double expect{};
		   double variance{};
		   //begining calculating the numerator of expectations and variances
		   for (size_t t = tau; t < m; t++){
			 for (size_t k = 0; k < (t-tau+1); k++) {
			   power1 = pow (phi, 2*k);
			   vecsum1.push_back(power1);
			 }
			 double sum1{};
			 for (auto& n : vecsum1)
			   sum1 += n;
			 numexp1 =  sum1 * pow(phi, (m-t)) * X[m];
			 vecsum1.clear();
			 for (size_t j = 0; j < (m-t); j++) {
			   power2 = pow (phi, 2*j);
			   vecsum2.push_back(power2);
			 }
			 double sum2{};
			 for (auto& n : vecsum2)
			   sum2 += n;
			 numexp2 = sum2 * pow(phi, (t-tau+1)) * X[tau-1];
			 vecsum2.clear();
			 variance = sigmasq*sum1*sum2 / den;
			 expect = (numexp1 + numexp2) / den;
			 expectations.push_back(expect);
			 variances.push_back(variance);
		   }
		   expectations.push_back(X[m]);
		   variances.push_back(0);
		   //endinging calculating the numerator of expectations and variances
		   missing.clear();
		 }
		 else {continue;}
		   }
		 }
		 //Create a dat files for expectations and variances
		 F_outExp(expectations);
		 F_outVar(variances);
		 cout << "case 3" << endl;
	}





	/****************case 4: No observation at the beginning None at the end ************/
	else {
		 // these are estimations of the first missing points. Step needed because
		 // I don't know the value of the first time
		 vector < size_t > missing{};
		 vector < double > expectations{};
		 vector < double > variances{};
		 size_t firstvalue{};
		 for ( size_t i = 0; i < N; i++ ){
		   if ( vect_obs_N[i] == 1 ) {
		 firstvalue = i;
		 break;
		   }
		 }
		 for ( size_t i = 0; i < firstvalue; i++ ){
		   double expect{};
		   double variance{};
		   double power{};
		   vector < double > vecsum{};
		   expect = pow(phi,(firstvalue-i)) * X[firstvalue];
		   power = pow(phi,2*(firstvalue-i-1));
		   expectations.push_back(expect);
		   vecsum.push_back(power);
		   double sum{};
		   for (auto& n : vecsum)
		 sum += n;
		   variance = sum*sigmasq;
		   variances.push_back(variance);
		 }
		 expectations.push_back(X[firstvalue]);
		 // end of the estimations of the first missing points.
		 size_t lastvalue{};
		 for ( int i = N-1; i >= 0; --i ){
		   if ( vect_obs_N[i] == 1 ) {
		 lastvalue = i+1;
		 break;
		   }
		 }
		 for ( size_t i = firstvalue+1; i < lastvalue; i++ ){
		   if ( (vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 1)) {
		 expectations.push_back(X[i]);
		 variances.push_back(0);
		   }
		   else if ((vect_obs_N[i] == vect_obs_N[i-1]) && (vect_obs_N[i] == 0)){}
		   else {
			 missing.push_back(i);
		 //these are all the calculations for variances and expectations
		 if (missing.size() == 2) {
		   size_t tau = missing[0];
		   size_t m = missing[1];
		   //beginning calculating the denomiator for variances and expectations
		   double den = F_den(m, tau, phi);
		   //ending calculating the denomiator for variances and expectations
		   double power1{};
		   double power2{};
		   vector < double > vecsum1{};
		   vector < double > vecsum2{};
		   double numexp1{};
		   double numexp2{};
		   double expect{};
		   double variance{};
		   //begining calculating the numerator of expectations and variances
		   for (size_t t = tau; t < m; t++){
			 for (size_t k = 0; k < (t-tau+1); k++) {
			   power1 = pow (phi, 2*k);
			   vecsum1.push_back(power1);
			 }
			 F_print_vector(vecsum1);
			 double sum1{};
			 for (auto& n : vecsum1)
			   sum1 += n;
			 numexp1 =  sum1 * pow(phi, (m-t)) * X[m];
			 vecsum1.clear();
			 for (size_t j = 0; j < (m-t); j++) {
			   power2 = pow (phi, 2*j);
			   vecsum2.push_back(power2);
			 }
			 F_print_vector(vecsum2);
			 double sum2{};
			 for (auto& n : vecsum2)
			   sum2 += n;
			 numexp2 = sum2 * pow(phi, (t-tau+1)) * X[tau-1];
			 vecsum2.clear();
			 variance = sigmasq*sum1*sum2 / den;
			 expect = (numexp1 + numexp2) / den;
			 expectations.push_back(expect);
			 variances.push_back(variance);
		   }
		   expectations.push_back(X[m]);
		   variances.push_back(0);
		   //endinging calculating the numerator of expectations and variances
		   missing.clear();
		 }
		 else {continue;}
		   }
		 }
		 // this is the estimation for the last missing points. It is needed because
		 // I don't know the value of the last time
		 for (size_t i = 0; i < N-lastvalue; i++) {
		 double expect{};
		 double variance{};
		 double power{};
		 vector < double > vecsum{};
		 expect = pow(phi,(i+1)) * X[lastvalue-1];
		 power = pow (phi,2*i);
		 expectations.push_back(expect);
		 vecsum.push_back(power);
		 double sum{};
		 for (auto& n : vecsum)
		   sum += n;
		 variance = sum*sigmasq;
		 variances.push_back(variance);
		 }
		 //Create a dat files for expectations and variances
		 F_outExp(expectations);
		 F_outVar(variances);
		 cout << "case 4" << endl;
	}

	return 0;

}



//*********** FUNCTIONS ***************

//Prints a matrix of unsigned size_t
void F_print_matrix ( vector < vector < size_t > > m ){
  for ( const vector < size_t > v : m ){
	for  ( size_t x : v ) cout << x << ' ';
	cout << endl;
  }
}
//Prints a matrix of signed doubles
void F_print_matrix ( vector < vector < double > > M ){
  for ( const vector < double > v : M ){
	for  ( double x : v ) cout << x << ' ';
	cout << endl;
  }
}
//Prints a vector of doubles
void F_print_vector ( vector < double > v ){
  for ( const double x : v ) cout << x << ' ';
  cout << endl;
  }

//Prints a vector of doubles
void F_print_vector(vector < size_t > v) {
	for (const size_t x : v) cout << x << ' ';
	cout << endl;
}

//Creates a dat file with the values of the Expectations
void F_outExp (vector <double> exp){
	ofstream outFile("./AR1_interp_exp_000.csv");
	outFile << endl;
	for (double n : exp) {
		outFile << n << endl;
	}
	outFile.close();
}

//Creates a dat file with the values of the Variances
void F_outVar(vector <double> var) {
	ofstream outFile("./AR1_interp_var_000.csv");
	outFile << endl;
	for (double n : var) {
		outFile << n << endl;
	}
	outFile.close();
}

//calculate the denominator for variance and expectation
double F_den(int m, int tau, double phi) {
	double den{};
	vector < double > vecsum{};
	double power{};
	for (size_t j = 0; j < (m - tau); j++) {
		power = pow(phi, 2 * j);
		vecsum.push_back(power);
	}
	for (auto& n : vecsum) den += n;
	return den;
}