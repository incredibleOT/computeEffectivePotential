#ifndef EFFECTIVEPOTENTIAL_H
#define EFFECTIVEPOTENTIAL_H

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

// Improved computation of the effective Potential
// U(v) = U_f(v) + 0.5*m_0^2*v^2 + lambda*v^4 + lambda_6*v^6 + BosDet
//       + v^2 * ( 6*lambda * (P_H + P_H) + lambda_6 * (45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ) +
//       + v^4 * lambda_6 *(15*P_H + 9*P_G)
//
// U_f(v) = -2*N_f/V * \sum_p { log( |nu^+(p) + y_t*v*(1 - 0.5*rho * nu^+(p))|^2 } - 2*N_f/V * \sum_p { log( |nu^+(p) + y_b*v*(1 - 0.5*rho * nu^+(p))|^2 }
// P_H = 1/V \sum_p { 1/( \hat{p}^2 + m_0^2 + 12*lambda*v^2 + 30 *lambda_6*v^4 ) } //NOTE without zero mode
// P_G = 1/V \sum_p { 1/( \hat{p}^2 + m_0^2 + 4*lambda*v^2 + 6 *lambda_6*v^4 ) } //NOTE without zero mode
//BosDet = 1/V \sum_p { 0.5 * log ( \hat{p}^2 + m_0^2 + 12*lambda*v^2 + 30 *lambda_6*v^4 ) + 1.5 * log ( \hat{p}^2 + m_0^2 + 4*lambda*v^2 + 6 *lambda_6*v^4 ) //NOTE without zero mode
//
// Approach to determin the Higgs boson mass is the following:
// set the lambdas and the desired cutoff -> vev
// Tune m_0^2 such that the 1st drivative of U wrt v is zero
//
//I'll try an iterative scheme where
// 1st: treelevel estimate is used to get first guess of m_0^2 ( m_0^2 =  - U_f'/v - 4*lambda*v^2 - 6*lambda_6*v^4
// iterate: insert m_0^2 into the propagatorsums and the determinant until a stable result is found
//
// An alternative approach would be: Set Goldstone mass to zero and determine higgs boson mass by the curvature of the potential at v
// determine iteratively m_0^2 until its stable NOTE: I don't know what to plug into the Determinant then
//
// NOTE: There is the problem that occurs, when m_0^2 gets large negative (such that \hat{p}^2 + m_0^2 + lambda*xxx get negative 

class effectivePotential
{
	private:
	int L0, L1, L2, L3;
	double cutoff_in_GeV;
	double y_t, y_b;
	double vev;
	double U_f, U_fp, U_fpp; //fermionic contributions at vev (are fixed for fixed, L, y and v)
	double propagatorSumWithZeroMass;//\sum{1/(\hat{p}^2} without (!!) zero momentum mode
	double lambda, lambda_6;
	
	//some constants
	double rho;
	double r; //hm, what is one somewhere else... with r=0.5, only r occurs in EV computation
	double vev_in_GeV;
	int N_f;
	
	//Propagator sums and derivatives wrt v and bosonic determinant with derivatives
	double actual_HPropSum, actual_HPropSumD1, actual_HPropSumD2;
	double actual_GPropSum, actual_GPropSumD1, actual_GPropSumD2;
	double actual_BosDet, actual_BosDetD1, actual_BosDetD2;
	
	//lattice_p^2 for bosons (4 * \sum_mu( sin^2( p_mu/2)))
	double smallest_pSquared; //for checking to not run into negative logarithms, stores the smallest non-zero values
	std::vector< double > pSquared;
	
	//sequence of m_0^2 and m_h^2 during the iteration
	std::vector< double > m0Squared, mHSquared;
	
	//maxCount for scans of potential
	static const int maxCount=10000;
	
	//================================
	void setConstants(); //sets (rho, r, N_f and vev_In_GeV);
	
	bool getFermionicContribution(); //assumes periodic BC so far
	
	void get_pSquared();
	
	
	std::complex< double > computeAnalyticalEigenvalue(double p0, double p1, double p2, double p3);
	
	// 	double computeBosonicDeterminant();
	
	void computeBosonicDeterminantAndPropagatorSums_WithDerivatives_combined(bool updateDeterminant);
	void computeBosonicDeterminantAndPropagatorSums_WithDerivatives();
	void computeBosonicDeterminantAndPropagatorSums_WithDerivatives_noDeterminantUpdate();//does not compute the log of the bosonic determinant, since it's not needed for the iteratrion
	
	//computes the bosonic determinant and its derivatives with m_g=0 amd mh^2 taken from the end of mHSquared;
// 	void computeBosonicDeterminantAndPropagatorSums_WithMassesByHand(bool updateDeterminant); NOTE nonsense

	public:
	effectivePotential(int l0, int l1, int l2, int l3); //can determine propagatorSum
	effectivePotential(int l0, int l1, int l2, int l3, double cutoff, double yt, double yb); //can determine fermionic contribution
	
	double computePropagatorSum(double massSquared);
	double computePropagatorSum();
		
	//sets values by hand
	void set_yukawa(double yt, double yb); //also (re)computes fermionic contribution, if cutoff (i.e. vev) is given;
	void set_cutoff(double cutoff); //also (re)computes fermionic contribution, if yukawa couplings are given;
	void set_lambda_6(double l_6);
	void set_lambda(double l);
	void set_rho(double new_rho);
	void set_r(double new_r);
	void set_N_f(int new_N_f);
	
	void resetConstants(double new_rho, double new_r, int new_N_f); //recomputes fermionic contr.
	
	void set_lambdaToStabilityBound(); //lambda \geq -(30/48 + 24 P_G)*lambda_6.  P_G=goldstone propagator sum (without zero mode)
	
	void set_MassSquared(double new_m0Squared, double new_mHSquared);//clears the mass vectors and sets the given values as first entries
	
	
	//just get the values
	double get_y_t();
	double get_y_b();
	double get_cutoff_in_GeV();
	double get_vev();
	double get_lambda_6();
	double get_lambda();
	double get_U_f();
	double get_U_fp();
	double get_U_fpp();
	std::vector< double > get_m0Squared_sequence();
	std::vector< double > get_mHSquared_sequence();
	double get_m0Squared();
	double get_mHSquared();
	double get_mH_in_GeV();
	
	
	//for the mesurement
	void initializeTreeLevel(); //sets m0Squared to its treelevel-value (zero loop order) and computes the zero order higgs mass and stores them in the vectors
	
	//uses the approach explained in the top 
	double iterateMass_withM0inPropSumsAndBosonicDet(); //returns the absolout(!) change in m0^2
	
	//uses the approach from philipps thesis (replace masses in the propagatorsums by zero or the curvature for goldstone and higgs) and no determinant
	double iterateMassDetermination_withMassInPropSumByHand(); //returns the sum of the absolut changes of m0Squared and mHSquared 
	
	//uses the last values of the masses and writes the potential and its first two derivatives to the stream
	bool plotPotentialToStream(double min, double max, double step, std::ostream &output, bool withInfo);
	private:
	bool fillVectorWithEigenvaluesAndPrefac( std::vector< std::pair< std::complex< double >, std::complex< double > > > &toBeFilled);
	public:
	
	//scans only the first derivative of the potential in the given range (cheaper, since no log ;) )
	//uses the last values of mHSquared and m0Squared
	bool scanPotential_firstDerivative_withMassInPropSumByHand(double min, double max, double step, std::map< double, double > &result_U_p);

};


#endif
