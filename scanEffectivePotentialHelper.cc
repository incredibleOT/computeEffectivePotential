#include "scanEffectivePotentialHelper.h"


void prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS )
{
	paraI["L0"]=-1;
	paraI["L1"]=-1;
	paraI["L2"]=-1;
	paraI["L3"]=-1;
	
	//directly set to default values
	paraI["N_f"]=1;
	paraD["rho"]=1.0;
	paraD["r"]=0.5;
	
	paraI["scan_cutoff_in_GeV"]=0;
	paraD["cutoff_in_GeV"]=-1.0;
	paraD["cutoff_in_GeV_min"]=-1.0;
	paraD["cutoff_in_GeV_max"]=-1.0;
	paraD["cutoff_in_GeV_step"]=-1.0;
	
	//if scan_y_t==0, y_t is choosen, if it is set, the range will be taken
	paraI["scan_y_t"]=0;
	paraD["y_t"]=-1.0;
	paraD["y_t_min"]=-1.0;
	paraD["y_t_max"]=-1.0;
	paraD["y_t_step"]=-1.0;
	
	//yRation=y_t/y_b
	paraI["scan_yRatio"]=0;
	paraD["yRatio"]=-1.0;
	paraD["yRatio_min"]=-1.0;
	paraD["yRatio_max"]=-1.0;
	paraD["yRatio_step"]=-1.0;
	
	
	paraI["scan_lambda_6"]=0;
	paraD["lambda_6"]=-1.0;
	paraD["lambda_6_min"]=-1.0;
	paraD["lambda_6_max"]=-1.0;
	paraD["lambda_6_step"]=-1.0;
	
	//determine_lambda==0 use scan and then like the others
	// ==1 stability bound (lambda>-lambda_6*(30/48 + 24*P_G))
	//==2 Lower lambda until there is a second minimum (consider the lowest lambda as the result with only one minimum in the potential in a given rang)
	// for the case 2: start at lambda_max and decrease by the abs of lambda_step
	//others may follow
	paraI["determine_lambda"]=-1;
	paraI["scan_lambda"]=0;
	paraD["lambda"]=-1.0;
	paraD["lambda_min"]=-1.0;
	paraD["lambda_max"]=-1.0;
	paraD["lambda_step"]=-1.0;
	
	//determin the iteration scheme
	//==0 use the one where the masses in the propagator sums are replaced by zero for goldstones and m_h for Higgs. Does not include bosonic determinant
	//==1 use (\hat{p}^2 + m_0^2 + 12*lambda*vev^2 + 30*lambda*vev^4) in propsums for higgs and respectively for golstoned  
	paraI["iteration_scheme"]=-1;
	
	paraS["OutputFile"] = "";
	
	//whether after determining the mass, the first derivative of the CEP schould be scanned
	//will test, whether a second extremum occurs and set the 
	paraI["scan_first_derivative"]=0;
	paraD["scan_derivative_min"]=-1.0;
	paraD["scan_derivative_max"]=-1.0;
	paraD["scan_derivative_step"]=-1.0;
	
	//whether the derivative should be printed to a file
	paraI["print_derivative_scan"]=0;
	//name will be: derivativeFileBody_ZZZZZ.txt
	//ZZZZZ will contain all quantities that were scanned
	//e.g.: _cut_x.x_; 
	paraS["derivativeFileBody"]="";
	
	
	
	
}


void prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet )
{
	prepareParameterMaps(paraD, paraI, paraS);
	for( std::map< std::string, double >::const_iterator iter=paraD.begin(); iter!=paraD.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
	for( std::map< std::string, int >::const_iterator iter=paraI.begin(); iter!=paraI.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
	for( std::map< std::string, std::string >::const_iterator iter=paraS.begin(); iter!=paraS.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
}


bool loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, const std::string &fileName )
{
	std::ifstream inputFile(fileName.c_str());
	std::string line,word;
	if(!(inputFile.good()))
	{
		std::cerr <<"Error opening inputfile " <<fileName <<std::endl;
		return false;
	}
	while(inputFile)
	{
		getline(inputFile, line);
		if(line.size()==0 || line[0] =='#' || line.find_first_not_of(' ') == std::string::npos){ continue; }
		std::istringstream strm(line);
		if(!(strm >> word)){ continue; }
		if( paraI.find(word) != paraI.end() ){ if(!(strm >> paraI[word])){ std::cerr <<"Error loading value of " <<word <<std::endl; std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; } }
		else if( paraD.find(word) != paraD.end() ){ if(!(strm >> paraD[word])){ std::cerr <<"Error loading value of " <<word <<std::endl; std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; } }
		else if( paraS.find(word) != paraS.end() ){ if(!(strm >> paraS[word])){ std::cerr <<"Error loading value of " <<word <<std::endl; std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; } }
		else{ std::cerr <<"Warning, there is no parameter called \"" <<word <<"\""<<std::endl; return false;}
	}
	inputFile.close(); inputFile.clear();
	return true;
}




bool loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, const std::string &fileName )
{
	std::ifstream inputFile(fileName.c_str());
	std::string line,word;
	if(!(inputFile.good()))
	{
		std::cerr <<"Error opening inputfile " <<fileName <<std::endl;
		return false;
	}
	while(inputFile)
	{
		getline(inputFile, line);
		if(line.size()==0 || line[0] =='#' || line.find_first_not_of(' ') == std::string::npos){ continue; }
		std::istringstream strm(line);
		if(!(strm >> word)){ continue; }
		if( paraI.find(word) != paraI.end() )
		{ 
			if(!(strm >> paraI[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else if( paraD.find(word) != paraD.end() )
		{ 
			if(!(strm >> paraD[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else if( paraS.find(word) != paraS.end() )
		{ 
			if(!(strm >> paraS[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else{ std::cerr <<"Warning, there is no parameter called \"" <<word <<"\""<<std::endl; return false;}
	}
	inputFile.close(); inputFile.clear();
	return true;
}




void streamParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::ostream &output, const std::string &prefix)
{
	for( std::map< std::string, double >::const_iterator iter=paraD.begin(); iter!=paraD.end(); ++iter )
	{
		if(!(prefix=="")){output <<prefix <<" ";}  
		output <<iter->first <<"     " <<iter->second <<std::endl;
	}
	for( std::map< std::string, int >::const_iterator iter=paraI.begin(); iter!=paraI.end(); ++iter )
	{
		if(!(prefix=="")){output <<prefix <<" ";}  
		output <<iter->first <<"     " <<iter->second <<std::endl;
	}
	for( std::map< std::string, std::string >::const_iterator iter=paraS.begin(); iter!=paraS.end(); ++iter )
	{
		if(!(prefix=="")){output <<prefix <<" ";}  
		output <<iter->first <<"     " <<iter->second <<std::endl;
	}
}



void streamSetParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, std::ostream &output, const std::string &prefix)
{
	for( std::map< std::string, double >::const_iterator iter=paraD.begin(); iter!=paraD.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
	for( std::map< std::string, int >::const_iterator iter=paraI.begin(); iter!=paraI.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
	for( std::map< std::string, std::string >::const_iterator iter=paraS.begin(); iter!=paraS.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
}




bool checkConsistencyOfParameters( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet )
{
	if( !( paraIsSet["L0"] && paraIsSet["L1"] && paraIsSet["L2"] && paraIsSet["L3"] ) )
	{
		std::cerr <<"Not all extends of the lattice are given" <<std::endl;
		return false;
	}
	if( paraI["L0"] <= 0 || paraI["L1"] <= 0 || paraI["L2"] <= 0 || paraI["L3"] <= 0 )
	{
		std::cerr <<"Not all extends of the lattice are positive" <<std::endl;
		return false;
	}
	//---------------------
	//cutoff_in_GeV
	if( paraIsSet["scan_cutoff_in_GeV"] && paraI["scan_cutoff_in_GeV"]!=0 )
	{
		if( !( paraIsSet["cutoff_in_GeV_min"] && paraIsSet["cutoff_in_GeV_max"] && paraIsSet["cutoff_in_GeV_step"] ) )
		{
			std::cerr <<"No scanrange in cutoff_in_GeV given" <<std::endl;
			return false;
		}
		if( paraD["cutoff_in_GeV_min"] >= paraD["cutoff_in_GeV_max"] || paraD["cutoff_in_GeV_step"] <= 0.0  )
		{
			std::cerr <<"Inconsistent scanning range in cutoff_in_GeV" <<std::endl;
			return false;
		}
		if( paraD["cutoff_in_GeV_min"] < 0.0 )
		{
			std::cerr <<"Negative cutoff_in_GeV is not allowed" <<std::endl;
			return false;
		}
	}
	//no scan
	if( !paraIsSet["scan_cutoff_in_GeV"] || paraI["scan_cutoff_in_GeV"]==0 )
	{
		if( !paraIsSet["cutoff_in_GeV"] || paraD["cutoff_in_GeV"] < 0.0 )
		{
			std::cerr <<"cutoff_in_GeV not set or negative" <<std::endl;
			return false;
		}
	}
	//----------------------
	//y_t - scan
	if( paraIsSet["scan_y_t"] && paraI["scan_y_t"]!=0 )
	{
		if( !( paraIsSet["y_t_min"] && paraIsSet["y_t_max"] && paraIsSet["y_t_step"] ) )
		{
			std::cerr <<"No scanrange in y_t given" <<std::endl;
			return false;
		}
		if( paraD["y_t_min"] >= paraD["y_t_max"] || paraD["y_t_step"] <= 0.0  )
		{
			std::cerr <<"Inconsistent scanning range in y_t" <<std::endl;
			return false;
		}
		if( paraD["y_t_min"] < 0.0 )
		{
			std::cerr <<"Negative y_t is not allowed" <<std::endl;
			return false;
		}
	}
	//no scan
	if( !paraIsSet["scan_y_t"] || paraI["scan_y_t"]==0 )
	{
		if( !paraIsSet["y_t"] || paraD["y_t"] < 0.0 )
		{
			std::cerr <<"y_t not set or negative" <<std::endl;
			return false;
		}
	}
	//----------------------
	//yRatio - scan
	if( paraIsSet["scan_yRatio"] && paraI["scan_yRatio"]!=0 )
	{
		if( !( paraIsSet["yRatio_min"] && paraIsSet["yRatio_max"] && paraIsSet["yRatio_step"] ) )
		{
			std::cerr <<"No scanrange in yRatio given" <<std::endl;
			return false;
		}
		if( paraD["yRatio_min"] >= paraD["yRatio_max"] || paraD["yRatio_step"] <= 0.0  )
		{
			std::cerr <<"Inconsistent scanning range in yRatio" <<std::endl;
			return false;
		}
		if( paraD["yRatio_min"] < 0.0 )
		{
			std::cerr <<"Negative yRatio is not allowed" <<std::endl;
			return false;
		}
	}
	//no scan
	if( !paraIsSet["scan_yRatio"] || paraI["scan_yRatio"]==0 )
	{
		if( !paraIsSet["yRatio"] || paraD["yRatio"] < 0.0 )
		{
			std::cerr <<"yRatio not set or negative" <<std::endl;
			return false;
		}
	}
	//----------------------
	//lambda_6 - scan
	if( paraIsSet["scan_lambda_6"] && paraI["scan_lambda_6"]!=0 )
	{
		if( !( paraIsSet["lambda_6_min"] && paraIsSet["lambda_6_max"] && paraIsSet["lambda_6_step"] ) )
		{
			std::cerr <<"No scanrange in lambda_6 given" <<std::endl;
			return false;
		}
		if( paraD["lambda_6_min"] >= paraD["lambda_6_max"] || paraD["lambda_6_step"] <= 0.0  )
		{
			std::cerr <<"Inconsistent scanning range in lambda_6" <<std::endl;
			return false;
		}
		if( paraD["lambda_6_min"] < 0.0 )
		{
			std::cerr <<"Negative lambda_6 is not allowed" <<std::endl;
			return false;
		}
	}
	//no scan
	if( !paraIsSet["scan_lambda_6"] || paraI["scan_lambda_6"]==0 )
	{
		if( !paraIsSet["lambda_6"] || paraD["lambda_6"] < 0.0 )
		{
			std::cerr <<"lambda_6 not set or negative" <<std::endl;
			return false;
		}
	}
	
	//------------------
	//lambda - determ
	if( ! paraIsSet["determine_lambda"] )
	{
		std::cerr <<"determine_lambda is not set" <<std::endl;
		return false;
	}
	if( paraI["determine_lambda"]==0)
	{
		//lambda - scan
		if( paraIsSet["scan_lambda"] && paraI["scan_lambda"]!=0 )
		{
			if( !( paraIsSet["lambda_min"] && paraIsSet["lambda_max"] && paraIsSet["lambda_step"] ) )
			{
				std::cerr <<"No scanrange in lambda given" <<std::endl;
				return false;
			}
			if( paraD["lambda_min"] >= paraD["lambda_max"] || paraD["lambda_step"] <= 0.0  )
			{
				std::cerr <<"Inconsistent scanning range in lambda" <<std::endl;
				return false;
			}
		}
		//no scan
		if( !paraIsSet["scan_lambda"] || paraI["scan_lambda"]==0 )
		{
			if( !paraIsSet["lambda"] )
			{
				std::cerr <<"lambda not set" <<std::endl;
				return false;
			}
		}
	}
	else if( paraI["determine_lambda"] == 1 )
	{
		//do nothing
	}
	else if( paraI["determine_lambda"] == 2 )
	{
		//check if lambda_max and lambda_step are set
		if( !paraIsSet["lambda_max"] || !paraIsSet["lambda_step"] || paraD["lambda_step"]==0.0 )
		{
			std::cerr <<"Set lambda_max and lambda_step for determine_lambda=2!" <<std::endl;
			return false;
		}
		if( !(paraIsSet["scan_first_derivative"] && paraI["scan_first_derivative"]!=0 ) )
		{
			std::cerr <<"Set scan_first_derivative for determine_lambda=2!" <<std::endl;
			return false;
		}
	}
	else
	{
		std::cerr <<"determine_lambda =" << paraI["determine_lambda"] <<" is not implemented" <<std::endl;
		return false;
	}
	//N_f
	if( paraIsSet["N_f"] && paraI["N_f"] <=0 )
	{
		std::cerr <<"N_f has to be positive if set" <<std::endl;
		return false;
	}
	//rho
	if( paraIsSet["rho"] && paraD["rho"] <=0 )
	{
		std::cerr <<"rho has to be positive if set" <<std::endl;
		return false;
	}
	//r
	if( paraIsSet["r"] && paraD["r"] <=0 )
	{
		std::cerr <<"r has to be positive if set" <<std::endl;
		return false;
	}
	//iteration_scheme
	if( !paraIsSet["iteration_scheme"] )
	{
		std::cerr <<"r has to be positive if set" <<std::endl;
		return false;
	}
	if( paraI["iteration_scheme"] < 0 || paraI["iteration_scheme"] > 1 )
	{
		std::cerr <<"iteration_scheme =" << paraI["iteration_scheme"] <<" is not implemented" <<std::endl;
		return false;	
	}
	//OutputFile
	if( !paraIsSet["OutputFile"] || paraS["OutputFile"]=="" )
	{
		std::cerr <<"OutputFile not set" <<std::endl;
		return false;
	}
	//scan_first_derivative
	if( paraIsSet["scan_first_derivative"] && paraI["scan_first_derivative"]!=0 )
	{
		if( !( paraIsSet["scan_derivative_min"] && paraIsSet["scan_derivative_max"] && paraIsSet["scan_derivative_step"] ) )
		{
			std::cerr <<"No scanrange in scan_derivative given" <<std::endl;
			return false;
		}
		if( paraD["scan_derivative_min"] >= paraD["scan_derivative_max"] || paraD["scan_derivative_step"] <= 0.0  )
		{
			std::cerr <<"Inconsistent scanning range in scan_derivative" <<std::endl;
			return false;
		}
		if( paraD["scan_derivative_min"] <= 0.0  )
		{
			std::cerr <<"Cannot scan derivative for non-positive vev" <<std::endl;
			return false;
		}
		//so far scan is only available for iteration_scheme=0
		if( paraI["iteration_scheme"]!=0 )
		{
			std::cerr <<"scan_derivative only implemented for iteration_scheme=0" <<std::endl;
			return false;
		}
	}
	if( paraIsSet["print_derivative_scan"] && paraI["print_derivative_scan"]!=0 )
	{
		if(!paraIsSet["derivativeFileBody"] || paraS["derivativeFileBody"]=="")
		{
			std::cerr <<"No derivativeFileBody given to print derivative scan" <<std::endl;
			return false;
		}
		if(!paraIsSet["scan_first_derivative"] || paraI["scan_first_derivative"]==0)
		{
			std::cerr <<"Output of derivative scan set, but no scan was set" <<std::endl;
			return false;
		}
	}
	return true;
}

void fillSetWithRange( const double min, const double max, const double step, std::set< double > &toFill)
{
	toFill.clear();
	
	int numberOfEntries=static_cast< int >( (max-min)/step + 1.5);
	
	for( int i=0; i<numberOfEntries; ++i)
	{
		toFill.insert( min + static_cast< double >(i) *step);
	}
	

}


bool printResultsVectorToStream(const std::vector< resultForOutput > &results, std::ostream &output)
{
	//increase precision, but store the old one
	std::streamsize oldPrec=output.precision();
	output.precision(12);
	for(std::vector< resultForOutput >::const_iterator iter=results.begin(); iter!=results.end(); ++iter)
	{
		if( !( output <<iter->cutoff_in_GeV <<" " <<iter->vev <<" " <<iter->y_t <<" " <<iter->y_b <<" " <<iter->lambda <<" "
		              <<iter->lambda_6 <<" " <<iter->m0Squared <<" " <<iter->mHSquared <<" " <<iter->mH_in_GeV <<" " <<iter->resultFlag <<std::endl ) )
		{
			std::cerr <<"Error during output of results" <<std::endl;
			return false;
		}
	}
	output.precision(oldPrec);
	return true;

}



int countNumberOfSignChanges(const std::map< double, double > &firstDerivative)
{
	//an exact zero also counts as signChange (even if its only a saddle point)
	int result=0;
	if(firstDerivative.empty() || firstDerivative.size()==1)
	{
		std::cerr <<"Error, empty map or only one element in countNumberOfSignchanges()"<<std::endl; return -1; 
	}
	std::map< double, double >::const_iterator iter=firstDerivative.begin();
	int signDummy(0);
	if( iter->second < 0.0 ){ signDummy=-1; }
	else if( iter->second > 0.0 ){ signDummy=1; }
	else{ signDummy=0; ++result; }
	++iter;
	while(iter!=firstDerivative.end())
	{
		if( iter->second < 0.0 )
		{ 
			if(signDummy==+1){ ++result; } 
			signDummy=-1; 
		}
		else if( iter->second > 0.0 )
		{ 
			if(signDummy==-1){ ++result; }
			signDummy=+1;
		}
		else
		{
			++result;
			signDummy=0;
		}
		++iter;
	}
	return result;
}

