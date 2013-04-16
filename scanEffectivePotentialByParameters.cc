//determines the higgs mass from the effective potential
//works with input file
//lattice extend is fixed
//scan possible in:
//y_t
//ratio y_t/y_b
//lambda_6
//lambda_0
//for lambda_0 there is are various possibilities to be set:
// - by hand
// - a stability criterion (\lambda > -\lambda_6*(30/48 + 24*P_G)
// - lower it, until there is a second extremum next to the minimum for the vev

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
	size_t numberOfParameters=parametersDouble.size()+parametersInt.size()+parametersString.size();
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
	
	
	cout <<"Parameters loaded:" <<endl;
	streamSetParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet, cout);
	
	cout <<"Check input parameters" <<endl;
	if( !checkConsistencyOfParameters(parametersDouble, parametersInt, parametersString, parametersIsSet) )
	{
		cerr <<"Error, given input parameters are inconsistent" <<endl;
		exit(EXIT_FAILURE);
	}
	else{ cout <<"passed" <<endl; }
	
	//vector to hold the results;
	std::vector< resultForOutput > results;
	
	
	//initiallize effPot
	bool useLowMemory(false);
	if( parametersInt["LowMemoryUsage"]!=0 ){ useLowMemory=true; }
	
		
	effectivePotential effPot(parametersInt["L0"], parametersInt["L1"], parametersInt["L2"], parametersInt["L3"], useLowMemory);
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
				//iterate lambda_6
				for(std::set< double >::const_iterator lambda_6=lambda_6_values.begin(); lambda_6!=lambda_6_values.end(); ++lambda_6)
				{
					cout <<"evaluating lambda_6=" <<*lambda_6 <<endl;
					effPot.set_lambda_6(*lambda_6);
					
					
					//set lambda
					std::set< double >lambda_values;
					
					resultForOutput dummy;
					
					//following needed for determine_lambda=2 or 3
					resultForOutput old_dummy;
					std::map< double, double > old_firstDerivativeOfPotential;
					std::map< double, double > old_secondDerivativeOfPotential;
					double lambda_step_actual=0.0;
					if( parametersInt["determine_lambda"]==0 && ( parametersInt["scan_lambda"]==0 ) )
					{ 
						lambda_values.insert( parametersDouble["lambda"] );
					}
					else if( parametersInt["determine_lambda"]==0 && parametersInt["scan_lambda"]!=0 )
					{ 
						fillSetWithRange( parametersDouble["lambda_min"], parametersDouble["lambda_max"], parametersDouble["lambda_step"], lambda_values );
					}
					else if( parametersInt["determine_lambda"]==1 )
					{
						lambda_values.insert( 0.0 );
					}
					else if( parametersInt["determine_lambda"]==2 || parametersInt["determine_lambda"]==3)
					{ 
						lambda_values.insert(parametersDouble["lambda_max"]);
						lambda_step_actual=std::abs(parametersDouble["lambda_step_max"]);
					}
					else
					{
						cerr <<"Error, unexpected case for setting lambda!" <<endl;
						exit(EXIT_FAILURE);
					}
					
					
					//iterate lambda
					std::set< double >::iterator lambda=lambda_values.begin();
					while(lambda != lambda_values.end())
					{
						if(parametersInt["determine_lambda"]==0)
						{
							effPot.set_lambda(*lambda);
							cout <<"evaluating lambda = " <<effPot.get_lambda() <<endl;
						}
						else if(parametersInt["determine_lambda"]==1)
						{
							effPot.set_lambdaToStabilityBound(); 
							cout <<"evaluating lambda = " <<effPot.get_lambda() <<" (from stability bound 1)" <<endl;
						}
						else if( parametersInt["determine_lambda"]==2 || parametersInt["determine_lambda"]==3 )
						{
							effPot.set_lambda( *lambda );
							cout <<"evaluating lambda = " <<effPot.get_lambda() <<endl;
						}
						
						
						int resultFlag(0);
						//now start the action
						effPot.initializeTreeLevel();
						cout <<"treelevel: m0Squared = " <<effPot.get_m0Squared() <<"   mHSquared = " <<effPot.get_mHSquared() <<endl;
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
									cout <<"convergence frozen" <<endl;
									resultFlag=1;
									break;
								}
								convergence_freeze_dummy=diff;
								resultFlag=1;
							}
						}
						
						
						dummy.cutoff_in_GeV = effPot.get_cutoff_in_GeV();
						dummy.vev           = effPot.get_vev();
						dummy.y_t           = effPot.get_y_t();
						dummy.y_b           = effPot.get_y_b();
						dummy.lambda        = effPot.get_lambda();
						dummy.lambda_6      = effPot.get_lambda_6() ;
						dummy.resultFlag    = resultFlag;
						if(resultFlag==2) //negative difference, hence error in determination
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
						
						//scan1stDerivative
						std::map< double, double > firstDerivativeOfPotential;
						int nExtreme=0;
						if( resultFlag!=2 &&  parametersInt["scan_first_derivative"] )
						{
							switch(parametersInt["iteration_scheme"])
							{
								case 0 :
											if( effPot.scanPotential_firstDerivative_withMassInPropSumByHand(parametersDouble["scan_derivative_min"], parametersDouble["scan_derivative_max"], parametersDouble["scan_derivative_step"], firstDerivativeOfPotential) )
											{
												nExtreme=countNumberOfSignChanges( firstDerivativeOfPotential );
												if (nExtreme!=1){ dummy.resultFlag=3; }
											}
											else{ resultFlag=2; }
								break;
								case 1 : break;
								default: cerr <<"Error, you should not be here" << endl; exit(EXIT_FAILURE);
							}
						}
						//scan2ndDerivative
						std::map< double, double > secondDerivativeOfPotential;
						bool negativeCurvatureAboveVev(true);
						if( resultFlag!=2 &&  parametersInt["scan_second_derivative"] )
						{
							switch(parametersInt["iteration_scheme"])
							{
								case 0 :
											if( effPot.scanPotential_secondDerivative_withMassInPropSumByHand(parametersDouble["scan_derivative_min"], parametersDouble["scan_derivative_max"], parametersDouble["scan_derivative_step"], secondDerivativeOfPotential) )
											{
												negativeCurvatureAboveVev=checkForNegativeCurvature(secondDerivativeOfPotential, effPot.get_vev());
											
												if(negativeCurvatureAboveVev){ dummy.resultFlag=3; }
											}
											else{ resultFlag=2; }
								break;
								case 1 : break;
								default: cerr <<"Error, you should not be here" << endl; exit(EXIT_FAILURE); 
							}
						}
						
						
						/*
						if( resultFlag!=2 &&  parametersInt["scan_first_derivative"] )
						{
							std::map< double, double > derivativeOfPotential;
							bool derivative_sucess(false);
							int nExtreme=0;
							switch(parametersInt["iteration_scheme"])
							{
								case 0 : derivative_sucess=effPot.scanPotential_firstDerivative_withMassInPropSumByHand(parametersDouble["scan_derivative_min"], parametersDouble["scan_derivative_max"], parametersDouble["scan_derivative_step"], derivativeOfPotential); nExtreme=countNumberOfSignChanges( derivativeOfPotential ); break;
								case 1 : derivative_sucess=false; break;
								default: derivative_sucess=false; 
							}
							
							if(derivative_sucess)
							{
								cout <<"First derivative shows " <<nExtreme <<" extrema"  <<endl;
								if(nExtreme!=1){ dummy.resultFlag=3; } 
								//NOTE holymagicshit
								//if there is no second extremum, store dummy in old_dummy and derivativeOfPotential
								//in old_derivativeOfPotential
								//further remove ++lambda_values.begin() and add lambda_values.begin()-abs(lambda_step)
								//set lambda to lambda_values.begin() and execute continue 
								//if there is a second extremum, lambda is too small. Copy old_dummy and
								//old_derivativeOfPotential and go on
								//I know, it's somehow bad style...
								if( parametersInt["determine_lambda"]==2 && dummy.mHSquared > 0.0 && nExtreme==1 )
								{
									old_dummy=dummy;
									old_derivativeOfPotential=derivativeOfPotential;
									lambda_values.erase( ++lambda_values.begin() );
									lambda_values.insert( *lambda_values.begin() - std::abs(parametersDouble["lambda_step"]) );
									cout <<"actual m0Sqared = " <<effPot.get_m0Squared() <<"     mHSquared = " <<effPot.get_mHSquared() <<"     mH = " <<effPot.get_mH_in_GeV() 
									<<"   valid. Decrease lambda." <<endl;
									lambda=lambda_values.begin();
									continue;
								}
								else if( parametersInt["determine_lambda"]==2 )
								{
									cout <<"Negative mass or multiple extrema, take last result!" <<endl;
									dummy=old_dummy;
									derivativeOfPotential=old_derivativeOfPotential;
									lambda=--(lambda_values.end());
								}
							}
							//output of derivative
							if(parametersInt["print_derivative_scan"]!=0)
							{
								std::ostringstream outputFileName;
								outputFileName<<parametersString["derivativeFileBody"];
								if(parametersInt["scan_cutoff_in_GeV"]!=0)
								{
									outputFileName<<"_cutoff_"<<dummy.cutoff_in_GeV;
								}
								if(parametersInt["scan_y_t"]!=0)
								{
									outputFileName<<"_yt_"<<dummy.y_t;
								}
								if(parametersInt["scan_yRatio"]!=0)
								{
									outputFileName<<"_yRatio_"<<dummy.y_t/dummy.y_b;
								}
								if(parametersInt["scan_lambda_6"]!=0)
								{
									outputFileName<<"_lambda6_"<<dummy.lambda_6;
								}
								if(parametersInt["determine_lambda"]==0 && parametersInt["scan_lambda"]!=0)
								{
									outputFileName<<"_lambda_"<<dummy.lambda;
								}
								outputFileName<<".txt";
								std::ofstream outputFile( outputFileName.str().c_str() );
								if(!outputFile.good())
								{
									cerr <<"Error opening output file" <<endl <<outputFileName.str() <<endl;
									exit(EXIT_FAILURE);
								}
								outputFile <<"# first derivative of the constrained effective Potential" <<endl;
								outputFile <<"# actual values used:" <<endl;
								outputFile <<"# L0: " <<parametersInt["L0"];
								outputFile <<"  L1: " <<parametersInt["L1"];
								outputFile <<"  L2: " <<parametersInt["L2"];
								outputFile <<"  L3: " <<parametersInt["L3"] <<endl;;
								outputFile <<"# Cutoff in GeV: " <<dummy.cutoff_in_GeV <<endl;
								outputFile <<"# y_t: " << dummy.y_t <<endl;
								outputFile <<"# y_b: " << dummy.y_b <<endl;
								outputFile <<"# lambda_6: " << dummy.lambda_6 <<endl;
								outputFile <<"# lambda: " << dummy.lambda <<endl;
								outputFile <<"# m0Squared: " << dummy.m0Squared <<endl;
								outputFile <<"# mHSquared: " << dummy.mHSquared <<endl;
								outputFile <<"# format is: vev   U_prime" <<endl;
								outputFile.precision(12);
								for(std::map< double, double >::const_iterator iter=derivativeOfPotential.begin(); iter!=derivativeOfPotential.end(); ++iter)
								{
									outputFile <<iter->first <<" " <<iter->second <<endl;
								}
								outputFile.close();
							}
						}*/
						
						if( parametersInt["determine_lambda"]==0 || parametersInt["determine_lambda"]==1 )
						{
							++lambda;
						}
						else if( parametersInt["determine_lambda"]==2 || parametersInt["determine_lambda"]==3)
						{
							bool goodLambda(false);
							if( parametersInt["determine_lambda"]==2 && nExtreme==1 ){ goodLambda=true; }
							else if(parametersInt["determine_lambda"]==3 && !negativeCurvatureAboveVev ){ goodLambda=true; }
							
							if(goodLambda)
							{
								//if lambda is good, copy result
								old_dummy=dummy;
								old_firstDerivativeOfPotential=firstDerivativeOfPotential;
								old_secondDerivativeOfPotential=secondDerivativeOfPotential;
								//set new lambda.
								//if lambda_step_actual is still lambda_step_max, take lambda_step_actual
								//if not, halve it again (as long, as lambda_step_actual is larger than lambda_step_min
								if(lambda_step_actual == std::abs(parametersDouble["lambda_step_max"]) )
								{
									//*lambda = old_dummy.lambda - lambda_step_actual; NOTE does not work on set!!
									lambda_values.clear();
									lambda_values.insert(old_dummy.lambda - lambda_step_actual);
									lambda=lambda_values.begin();
									continue; //use next lambda but skip the result storage
								}
								else if(lambda_step_actual >= parametersDouble["lambda_step"])
								{
									lambda_step_actual/=2.0;
									//*lambda = old_dummy.lambda - lambda_step_actual; NOTE does not work on set!!
									lambda_values.clear();
									lambda_values.insert(old_dummy.lambda - lambda_step_actual);
									lambda=lambda_values.begin();
									continue; //use next lambda but skip the result storage
								}
								else
								{
									// smallest stepsize reached, print result, i.e. do nothing more 
									// right dummy is used anyway
									++lambda; //to leave loop
								}
							}
							else if(old_dummy.resultFlag != -1) //resultFlag=-1 means no old results available
							{
								//old_dummy.resultFlag==-1 means no old results available, i.e. already the first lambda resulted 
								//in a bad potential
								//
								//if lambda is not good, but stepsize not yet the smallest decrease lambda_step_actual
								//take last working lambda  half of the stepsize 
								if(lambda_step_actual >= parametersDouble["lambda_step"])
								{
									lambda_step_actual/=2.0;
									//*lambda = old_dummy.lambda - lambda_step_actual; NOTE does not work on set!!
									lambda_values.clear();
									lambda_values.insert(old_dummy.lambda - lambda_step_actual);
									lambda=lambda_values.begin();
									continue;  //use next lambda but skip the result storage
								}
								else
								{
									//take old result
									dummy=old_dummy;
									firstDerivativeOfPotential=old_firstDerivativeOfPotential;
									secondDerivativeOfPotential=old_secondDerivativeOfPotential;
									++lambda; //to leave loop
								}
							}
							else
							{
								//already the first try resulted in a bad potential, store actual result
								// resultFlag will indicate, that there is a problem (it should be 3 or 2)
								++lambda; //to leave loop
							}
						}
						
						
						//print derivatives to file if desired
						if(parametersInt["print_derivative_scan"]!=0)
						{
							for(int whichDerivative=1; whichDerivative<=2; ++whichDerivative)
							{
								if(whichDerivative==1 && parametersInt["scan_first_derivative"]==0)
								{
									continue;
								}
								else if(whichDerivative==2 && parametersInt["scan_second_derivative"]==0)
								{
									continue;
								}
								std::ostringstream outputFileName;
								outputFileName<<parametersString["derivativeFileBody"];
								if(parametersInt["scan_cutoff_in_GeV"]!=0)
								{
									outputFileName<<"_cutoff_"<<dummy.cutoff_in_GeV;
								}
								if(parametersInt["scan_y_t"]!=0)
								{
									outputFileName<<"_yt_"<<dummy.y_t;
								}
								if(parametersInt["scan_yRatio"]!=0)
								{
									outputFileName<<"_yRatio_"<<dummy.y_t/dummy.y_b;
								}
								if(parametersInt["scan_lambda_6"]!=0)
								{
									outputFileName<<"_lambda6_"<<dummy.lambda_6;
								}
								if(parametersInt["determine_lambda"]==0 && parametersInt["scan_lambda"]!=0)
								{
									outputFileName<<"_lambda_"<<dummy.lambda;
								}
								if(whichDerivative==1)
								{
									outputFileName<<"_1stDerivative";
								}
								else if(whichDerivative==2)
								{
									outputFileName<<"_2ndDerivative";
								}
								outputFileName<<".txt";
								std::ofstream outputFile( outputFileName.str().c_str() );
								if(!outputFile.good())
								{
									cerr <<"Error opening output file" <<endl <<outputFileName.str() <<endl;
									exit(EXIT_FAILURE);
								}
								if(whichDerivative==1)
								{
									outputFile <<"# first derivative of the constrained effective Potential" <<endl;
								}
								else if(whichDerivative==2)
								{
									outputFile <<"# second derivative of the constrained effective Potential" <<endl;
								}
								outputFile <<"# actual values used:" <<endl;
								outputFile <<"# L0: " <<parametersInt["L0"];
								outputFile <<"  L1: " <<parametersInt["L1"];
								outputFile <<"  L2: " <<parametersInt["L2"];
								outputFile <<"  L3: " <<parametersInt["L3"] <<endl;;
								outputFile <<"# Cutoff in GeV: " <<dummy.cutoff_in_GeV <<endl;
								outputFile <<"# y_t: " << dummy.y_t <<endl;
								outputFile <<"# y_b: " << dummy.y_b <<endl;
								outputFile <<"# lambda_6: " << dummy.lambda_6 <<endl;
								outputFile <<"# lambda: " << dummy.lambda <<endl;
								outputFile <<"# m0Squared: " << dummy.m0Squared <<endl;
								outputFile <<"# mHSquared: " << dummy.mHSquared <<endl;
								outputFile <<"# format is: vev   derivative" <<endl;
								outputFile.precision(12);
								if(whichDerivative==1)
								{
									for(std::map< double, double >::const_iterator iter=firstDerivativeOfPotential.begin(); iter!=firstDerivativeOfPotential.end(); ++iter)
									{
										outputFile <<iter->first <<" " <<iter->second <<endl;
									}
								}
								else if(whichDerivative==2)
								{
									for(std::map< double, double >::const_iterator iter=secondDerivativeOfPotential.begin(); iter!=secondDerivativeOfPotential.end(); ++iter)
									{
										outputFile <<iter->first <<" " <<iter->second <<endl;
									}
								}
								
								outputFile.close();
							}
						}
						results.push_back(dummy);
						cout <<"Result: m0Sqared = " <<dummy.m0Squared <<"     mHSquared = " <<dummy.mHSquared 
						     <<"     mH = " <<sqrt(dummy.mHSquared)*dummy.cutoff_in_GeV <<endl;
						
						
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
		outputFile <<"#              3=scan of potential showed not only one extremum" <<endl;
		if(! printResultsVectorToStream( results, outputFile ) )
		{
			cerr <<"Error, during output to" <<endl <<parametersString["OutputFile"] <<endl;
		}
		outputFile.close();
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
	
