//determines the higgs mass from the effective potential
//works with input file
//lattice extend is fixed
//scan possible in:
//y_t
//ratio y_t/y_b
//lambda_6
//lambda_0
//for lambda_0 there is also the possibility to set it by a stability criterion (choice of multiple criterions may be possible at some point)
// --------


#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
// #include <ofstream>
#include <ostream>
#include <set>
#include <sstream>
#include <string>

#include "effectivePotential.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
	if(narg!=11)
	{
		cerr <<"Error, start program with:" <<endl;
		cerr <<arg[0] <<"   L_s   L_t   cutoff   y_t   lambda_6   lambda   v_min   v_max   v_step   outputFile" <<endl;
		cerr <<"quark masses are degenerate" <<endl;
		exit(EXIT_FAILURE);
	}
	int L_s              = atoi(arg[1]);
	int L_t              = atoi(arg[2]);
	double cutoff_in_GeV = atof(arg[3]);
	double y_t           = atof(arg[4]);
	double lambda_6      = atof(arg[5]);
	double lambda        = atof(arg[6]);
	double v_min         = atof(arg[7]);
	double v_max         = atof(arg[8]);
	double v_step        = atof(arg[9]);
	std::string outputFileName = (arg[10]);
	
	cout <<"L_s = " <<L_s <<endl;
	cout <<"L_t = " <<L_t <<endl;
	cout <<"cutoff_in_GeV = " <<cutoff_in_GeV <<endl;
	cout <<"y_t = " <<y_t <<endl;
	cout <<"lambda_6 = " <<lambda_6 <<endl;
	cout <<"lambda = " <<lambda <<endl;
	cout <<"v_min = " <<v_min <<endl;
	cout <<"v_max = " <<v_max <<endl;
	cout <<"v_step = " <<v_step <<endl;
	cout <<"outputFileName = " <<outputFileName <<endl;
	
	effectivePotential effPot(L_s, L_s, L_s, L_t, cutoff_in_GeV, y_t, y_t);
	effPot.set_lambda_6(lambda_6);
	effPot.set_lambda(lambda);
	
	effPot.initializeTreeLevel();
	double convergence_freeze_dummy(0.0);
	bool errorFlag(false);
	for(int i=1; i<=500; ++i)
	{
		double diff=-1.0;
		diff=effPot.iterateMassDetermination_withMassInPropSumByHand();
		if(diff <0.0)
		{
			cerr <<"Error, iteration in mass determination results in negative absolute mass difference" <<endl;
			errorFlag=true;
			break;
		}
		if(diff < 1e-16)
		{
			cout <<"Massdetermination converged after " <<i <<" iteration. Last difference: " <<diff <<endl;
			errorFlag=false;
			break;
		}
		if( (i%25)==0 )
		{
			cout <<"Massdetermination did not converged after " <<i <<" iterations. Last difference: " <<diff <<endl;
			if(diff==convergence_freeze_dummy)
			{
				cout <<"convergence frozen" <<endl;
				break;
			}
			convergence_freeze_dummy=diff;
		}
	}
	if(!errorFlag)
	{
		std::ofstream outputFile(outputFileName.c_str());
		if(!(outputFile.good()))
		{
			cerr <<"Error opening" <<endl <<outputFileName <<endl;
			exit(EXIT_FAILURE);
		}
		effPot.plotPotentialToStream_withMassInPropSumByHand(v_min,v_max,v_step,outputFile,true);
		if(!(outputFile.good()))
		{
			cerr <<"Error during output to" <<endl <<outputFileName <<endl;
		}
		outputFile.close();
	}
}
	
	
