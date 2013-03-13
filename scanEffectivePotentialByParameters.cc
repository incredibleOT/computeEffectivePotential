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
#include "scanEffectivePotentialHelper.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
	std::map< std::string, double > parametersDouble;
	std::map< std::string, int > parametersInt;
	std::map< std::string, std::string > parametersString;
	std::map< std::string, bool > parametersIsSet;
	
	prepareParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet);
	if(narg!=2)
	{
		cerr <<"Error, start program with:" <<endl;
		cerr <<arg[0] <<"   inputFile" <<endl;
		streamParameterMaps(parametersDouble, parametersInt, parametersString, cerr);
		exit(EXIT_FAILURE);
	}
	if( !loadParameterMapsFromFile(parametersDouble, parametersInt, parametersString, parametersIsSet, arg[1]) )
	{
		cerr <<"Error loading input file" <<endl <<arg[1] <<endl;
		exit(EXIT_FAILURE);
	}
	size_t numberOfParameters=parametersDouble.size()+parametersInt.size()+parametersString.size();
	
	cout <<"Parameters loaded:" <<endl;
	streamSetParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet, cout);
	
	cout <<"Check input parameters" <<endl;
	if( !checkConsistencyOfParameters(parametersDouble, parametersInt, parametersString, parametersIsSet) )
	{
		cerr <<"Error, given input parameters are inconsitsnet" <<endl;
		exit(EXIT_FAILURE);
	}
	else{ cout <<"passed" <<endl; }
	
	//vector to hold the results;
	std::vector< resultForOutput > results;
	
	
	//initiallize effPot
	effectivePotential effPot(parametersInt["L0"], parametersInt["L1"], parametersInt["L2"], parametersInt["L3"]);
// 	effPot.set_cutoff( parametersDouble["cutoff_in_GeV"] );
	if(parametersIsSet["N_f"]){ effPot.set_N_f(parametersInt["N_f"]); }
	if(parametersIsSet["rho"]){ effPot.set_rho(parametersDouble["rho"]); }
	if(parametersIsSet["r"]){ effPot.set_r(parametersDouble["r"]); }
	
	//set cutoff_in_GeV 
	std::set< double > cutoff_in_GeV_values;if( !parametersIsSet["scan_cutoff_in_GeV"] || parametersInt["scan_cutoff_in_GeV"]==0 ){ cutoff_in_GeV_values.insert( parametersDouble["cutoff_in_GeV"] ); }
	else{ fillSetWithRange( parametersDouble["cutoff_in_GeV_min"], parametersDouble["cutoff_in_GeV_max"], parametersDouble["cutoff_in_GeV_step"], cutoff_in_GeV_values); }
	//iterate over cutoff_in_GeV
	for(std::set< double >::const_iterator cutoff_in_GeV=cutoff_in_GeV_values.begin(); cutoff_in_GeV!=cutoff_in_GeV_values.end(); ++cutoff_in_GeV)
	{
		cout <<"evaluating cutoff_in_GeV=" <<*cutoff_in_GeV <<endl;
		effPot.set_cutoff( *cutoff_in_GeV );
		//set y_t
		std::set< double > y_t_values;
		if( !parametersIsSet["scan_y_t"] || parametersInt["scan_y_t"]==0 ){ y_t_values.insert( parametersDouble["y_t"] ); }
		else{ fillSetWithRange( parametersDouble["y_t_min"], parametersDouble["y_t_max"], parametersDouble["y_t_step"], y_t_values); }
		//iterate over y_t
		for(std::set< double >::const_iterator y_t=y_t_values.begin(); y_t!=y_t_values.end(); ++y_t)
		{
			cout <<"evaluating y_t=" <<*y_t <<endl;
			
			//set yRatio
			std::set< double > yRatio_values;
			if( !parametersIsSet["scan_yRatio"] || parametersInt["scan_yRatio"]==0 ){ yRatio_values.insert( parametersDouble["yRatio"] ); }
			else{ fillSetWithRange( parametersDouble["yRatio_min"], parametersDouble["yRatio_max"], parametersDouble["yRatio_step"], yRatio_values ); }
			//iterate yRatio
			for(std::set< double >::const_iterator yRatio=yRatio_values.begin(); yRatio!=yRatio_values.end(); ++yRatio)
			{
				double y_b=(*y_t)/(*yRatio);
				cout <<"evaluating yRatio=" <<*yRatio <<" (y_b=" <<y_b <<")" <<endl;
				//set yukawa couplings
				effPot.set_yukawa(*y_t,y_b);
				
				//set lambda_6
				
				std::set< double > lambda_6_values;
				if( !parametersIsSet["scan_lambda_6"] || parametersInt["scan_lambda_6"]==0 ){ lambda_6_values.insert( parametersDouble["lambda_6"] ); }
				else{ fillSetWithRange( parametersDouble["lambda_6_min"], parametersDouble["lambda_6_max"], parametersDouble["lambda_6_step"], lambda_6_values ); }
				//cout <<"lambda_6_values.size()=" <<lambda_6_values.size() <<endl;
				//iterate lambda_6
				for(std::set< double >::const_iterator lambda_6=lambda_6_values.begin(); lambda_6!=lambda_6_values.end(); ++lambda_6)
				{
					cout <<"evaluating lambda_6=" <<*lambda_6 <<endl;
					effPot.set_lambda_6(*lambda_6);
					
					//set lambda
					bool useStabilityCriterion_1(false); // lambda=-lambda_6*(30/48 + 24* P_G)
					std::set< double >lambda_values;
					if( parametersInt["determine_lambda"]==0 && ( !parametersIsSet["scan_lambda"] || parametersInt["scan_lambda"]==0 ) ){ lambda_values.insert( parametersDouble["lambda"] ); }
					else if( parametersInt["determine_lambda"]==1 ){ useStabilityCriterion_1=true; lambda_values.insert( 0.0 ); }
					else{ fillSetWithRange( parametersDouble["lambda_min"], parametersDouble["lambda_max"], parametersDouble["lambda_step"], lambda_values ); }
					//iterate lambda
					for(std::set< double >::const_iterator lambda=lambda_values.begin(); lambda!=lambda_values.end(); ++lambda)
					{
						if(useStabilityCriterion_1){ effPot.set_lambdaToStabilityBound(); cout <<"evaluating lambda = " <<effPot.get_lambda() <<" (from stability bound 1)" <<endl;}
						else{ effPot.set_lambda(*lambda); cout <<"evaluating lambda = " <<effPot.get_lambda() <<endl; }
						
						int resultFlag(0);
						//now start the action
						effPot.initializeTreeLevel();
						double convergence_freeze_dummy(0.0);
						for(int i=1; i<=500; ++i)
						{
							double diff=-1.0;
							switch(parametersInt["iteration_scheme"])
							{
								case 0 : diff=effPot.iterateMassDetermination_withMassInPropSumByHand(); break;
								case 1 : diff=effPot.iterateMass_withM0inPropSumsAndBosonicDet(); break;
								default: diff=-1.0;
							}
							cout <<"iteration nr: " <<i <<", diff: " <<diff <<endl;
							if(diff <0.0)
							{
								cerr <<"Error, iteration in mass determination results in negative absolute mass difference" <<endl;
								resultFlag=2;
								break;
							}
							if(diff < 1e-16)
							{
								cout <<"Massdetermination converged after " <<i <<" iteration. Last difference: " <<diff <<endl;
								resultFlag=0;
								break;
							}
							if( (i%25)==0 )
							{
								cout <<"Massdetermination did not converged after " <<i <<" iterations. Last difference: " <<diff <<endl;
								if(diff==convergence_freeze_dummy)
								{
									cout <<"convergence freezed" <<endl;
									resultFlag=1;
									break;
								}
								convergence_freeze_dummy=diff;
								resultFlag=1;
							}
						}
						resultForOutput dummy;
						dummy.cutoff_in_GeV = effPot.get_cutoff_in_GeV();
						dummy.vev           = effPot.get_vev();
						dummy.y_t           = effPot.get_y_t();
						dummy.y_b           = effPot.get_y_b();
						dummy.lambda        = effPot.get_lambda();
						dummy.lambda_6      = effPot.get_lambda_6() ;
						dummy.resultFlag    = resultFlag;
						if(resultFlag==2)
						{
							dummy.m0Squared     = 0.0;
							dummy.mHSquared     = 0.0;
							dummy.mH_in_GeV     = 0.0;
						}
						else
						{
							dummy.m0Squared     = effPot.get_m0Squared();
							dummy.mHSquared     = effPot.get_mHSquared();
							dummy.mH_in_GeV     = effPot.get_mH_in_GeV();
						}
						
						results.push_back(dummy);
						cout <<"Result: m0Sqared = " <<effPot.get_m0Squared() <<"     mHSquared = " <<effPot.get_mHSquared() <<"     mH = " <<effPot.get_mH_in_GeV() <<endl;
					}
				}
			}
		}
	}
	
	//output
	{
		std::ostream &output=cout;
		if(! printResultsVectorToStream( results, output ) )
		{
			cerr <<"Error, during output" <<endl;
		}
		std::ofstream outputFile( parametersString["OutputFile"].c_str() );
		if(!outputFile.good())
		{
			cerr <<"Error opening output file" <<endl <<parametersString["OutputFile"] <<endl;
			exit(EXIT_FAILURE);
		}
		//info
// 		output=outputFile;
		outputFile <<"# Higgs bososn mass bounds from effective potential" <<endl;
		outputFile <<"# Set parameters:" <<endl;
		streamSetParameterMaps( parametersDouble, parametersInt, parametersString, parametersIsSet, outputFile, "#");
		outputFile <<"# Output format is:" <<endl;
		outputFile <<"# cutoff_in_GeV   vev   y_t   y_b   lambda   lambda_6   m0Squared   mHSquared   mH_in_GeV   resultFlag" <<endl;
		outputFile <<"# resultFlag:  0=no error   1=determination did not converge   2=Error during iteration" <<endl;
		if(! printResultsVectorToStream( results, outputFile ) )
		{
			cerr <<"Error, during output to" <<endl <<parametersString["OutputFile"] <<endl;
		}
	}
	
	
	
	if(parametersDouble.size()+parametersInt.size()+parametersString.size() != numberOfParameters || parametersIsSet.size()!=numberOfParameters )
	{
		cerr <<"Error, number of parameters changed!" <<endl;
		cerr <<"parametersDouble.size()=" <<parametersDouble.size() <<endl;
		cerr <<"parametersInt.size()=" <<parametersInt.size() <<endl;
		cerr <<"parametersString.size()=" <<parametersString.size() <<endl;
		cerr <<"parametersIsSet.size()=" <<parametersIsSet.size() <<endl;
		cerr <<"Iterating parametersIsSet:" <<endl;
		for(std::map< std::string, bool >::const_iterator iter=parametersIsSet.begin(); iter!=parametersIsSet.end(); ++iter)
		{
			cerr <<iter->first <<"   " <<iter->second <<endl;
		}
		exit(EXIT_FAILURE);
	}
	
}	
	
