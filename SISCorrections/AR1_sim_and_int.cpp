#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "SIS_AR1.h"

using namespace std;

void print_matrix(vector < vector < size_t > > m);
void print_matrix(vector < vector < double > > M);
void print_vector(vector < double > v);



int main() {

	random_device rd;
	mt19937 generator(rd());

	sis();


	/* This is the code to find the analytic estimates of expectation and variance
	   of the missing values in the AR(1) model. These estimates will then be used as
	   gold standard to evaluate our method.*/
	   
	   for ( size_t j = 0; j < N; j++) {cout << vect_obs_N[j] << endl;}
	   //case 1: I have one observation at the beginning and one at the end
	   if ( (vect_obs_N[0] == 1) && (vect_obs_N[N-1] == 1) ) {
		 vector < size_t > missing{};
		 vector < double > expectations{X[0]};
		 vector < double > variances{0};
		 for ( size_t i = 1; i < N; i++ ){
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
		   double den{};
		   vector < double > vecsum{};
		   double power{};
		   for ( size_t j = 0; j < (m-tau+1); j++) {
			 power = pow(phi, 2*j);
			 vecsum.push_back(power);
		   }
		   for (auto& n : vecsum)
			 den += n;
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
		 //Create a dat file with the values of the Expectations
		 ofstream outFile7( "./AR1_interp_exp.dat" );
		 outFile7 << endl;
		 for ( double n : expectations ){
		   outFile7 << n << endl;
		 }
		 outFile7.close();
		 //Create a dat file with the values of the Variances
		 ofstream outFile8( "./AR1_interp_var.dat" );
		 outFile8 << endl;
		 for ( double n : variances ){
		   outFile8 << n << endl;
		 }
		 outFile8.close();
	   }

	   //case 2: I have one observation at the beginning and none at the end
	   else if ( (vect_obs_N[0] == 1) && (vect_obs_N[N-1] == 0) ) {
		 vector < size_t > missing{};
		 vector < double > expectations{X[0]};
		 vector < double > variances{0};
		 size_t lastvalue{};
		 for ( int i = N-1; i >= 0; --i ){
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
			 missing.push_back(i);
		 //these are all the calculations for variances and expectations
		 if (missing.size() == 2) {
		   size_t tau = missing[0];
		   size_t m = missing[1];
		   //beginning calculating the denomiator for variances and expectations
		   double den{};
		   vector < double > vecsum{};
		   double power{};
		   for ( size_t j = 0; j < (m-tau+1); j++) {
			 power = pow(phi, 2*j);
			 vecsum.push_back(power);
		   }
		   for (auto& n : vecsum)
			 den += n;
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
		 // these are the estimations of the last missing points. Step needed because
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
		 // end of the estimations of the last missing points.
		 //Create a dat file with the values of the Expectations
		 ofstream outFile7( "./AR1_interp_exp.dat" );
		 outFile7 << endl;
		 for ( double n : expectations ){
		   outFile7 << n << endl;
		 }
		 outFile7.close();
		 //Create a dat file with the values of the Variances
		 ofstream outFile8( "./AR1_interp_var.dat" );
		 outFile8 << endl;
		 for ( double n : variances ){
		   outFile8 << n << endl;
		 }
		 outFile8.close();
	   }

	   //case 3: I have no observation at the beginning but one at the end
	   else if ( (vect_obs_N[0] == 0) && (vect_obs_N[N-1] == 1) ) {
		 // these are estimations of the first missing point. Step needed because
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
		   double den{};
		   vector < double > vecsum{};
		   double power{};
		   for ( size_t j = 0; j < (m-tau+1); j++) {
			 power = pow(phi, 2*j);
			 vecsum.push_back(power);
		   }
		   for (auto& n : vecsum)
			 den += n;
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
		 //Create a dat file with the values of the Expectations
		 ofstream outFile7( "./AR1_interp_exp.dat" );
		 outFile7 << endl;
		 for ( double n : expectations ){
		   outFile7 << n << endl;
		 }
		 outFile7.close();

		 //Create a dat file with the values of the Variances
		 ofstream outFile8( "./AR1_interp_var.dat" );
		 outFile8 << endl;
		 for ( double n : variances ){
		   outFile8 << n << endl;
		 }
		 outFile8.close();
	   }
	   //case 4: I have no observation at the beginning or at the end
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
		 cout << "first value " << firstvalue << endl;
		 expectations.push_back(X[firstvalue]);
		 // end of the estimations of the first missing points.
		 size_t lastvalue{};
		 for ( int i = N-1; i >= 0; --i ){
		   if ( vect_obs_N[i] == 1 ) {
		 lastvalue = i+1;
		 break;
		   }
		 }
		 cout << "last value " << lastvalue << endl;
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
		   cout << "tau" << tau << "m" << m << endl;
		   //beginning calculating the denomiator for variances and expectations
		   double den{};
		   vector < double > vecsum{};
		   double power{};
		   for ( size_t j = 0; j < (m-tau+1); j++) {
			 power = pow(phi, 2*j);
			 vecsum.push_back(power);
		   }
		   for (auto& n : vecsum)
			 den += n;
		   cout << " denominator " << den << endl;
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
			 print_vector(vecsum1);
			 double sum1{};
			 for (auto& n : vecsum1)
			   sum1 += n;
			 numexp1 =  sum1 * pow(phi, (m-t)) * X[m];
			 vecsum1.clear();
			 cout << "numexp1 " << numexp1 << endl;
			 cout << "sum1 " << sum1 << endl;
			 for (size_t j = 0; j < (m-t); j++) {
			   power2 = pow (phi, 2*j);
			   vecsum2.push_back(power2);
			 }
			 print_vector(vecsum2);
			 double sum2{};
			 for (auto& n : vecsum2)
			   sum2 += n;
			 numexp2 = sum2 * pow(phi, (t-tau+1)) * X[tau-1];
			 vecsum2.clear();
			 variance = sigmasq*sum1*sum2 / den;
			 cout << "numexp2 " << numexp2 << endl;
			 cout << "sum2 " << sum2 << endl;
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
		 print_vector(expectations);
		 print_vector(variances);
		 //Create a dat file with the values of the Expectations
		 ofstream outFile7( "./AR1_interp_exp.dat" );
		 outFile7 << endl;
		 for ( double n : expectations ){
		   outFile7 << n << endl;
		 }
		 outFile7.close();

		 //Create a dat file with the values of the Variances
		 ofstream outFile8( "./AR1_interp_var.dat" );
		 outFile8 << endl;
		 for ( double n : variances ){
		   outFile8 << n << endl;
		 }
		 outFile8.close();

	   }

	return 0;

}



//functions definitions
/*
//this function prints a matrix of unsigned size_t
void print_matrix ( vector < vector < size_t > > m ){
  for ( const vector < size_t > v : m ){
	for  ( size_t x : v ) cout << x << ' ';
	cout << endl;
  }
}
//this function prints a matrix of signed doubles
void print_matrix ( vector < vector < double > > M ){
  for ( const vector < double > v : M ){
	for  ( double x : v ) cout << x << ' ';
	cout << endl;
  }
}
//this function prints a vector of doubles
void print_vector ( vector < double > v ){
  for ( const double x : v ) cout << x << ' ';
  cout << endl;
  }*/