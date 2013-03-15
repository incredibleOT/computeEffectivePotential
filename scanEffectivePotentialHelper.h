#ifndef SCANEFFECTIVEPOTENTIALHERLPER_H
#define SCANEFFECTIVEPOTENTIALHERLPER_H


#include <istream>
#include <iostream>
#include <fstream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

void prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS );

void prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );

bool loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, const std::string &fileName );

bool loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, const std::string &fileName );

void streamParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::ostream &output, const std::string &prefix = "");

void streamSetParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, std::ostream &output, const std::string &prefix = "");

bool checkConsistencyOfParameters( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );

void fillSetWithRange( const double min, const double max, const double step, std::set< double > &toFill );

struct resultForOutput
{
	double cutoff_in_GeV;
	double vev;
	double y_t;
	double y_b;
	double lambda;
	double lambda_6;
	double m0Squared;
	double mHSquared;
	double mH_in_GeV;
	int resultFlag; //0:no error; 1:no convergence(frosen or maxcount reached), 2:error during iteration, 3:not only 1 extremum; 
	resultForOutput(): cutoff_in_GeV(0.0), vev(0.0), y_t(0.0), y_b(0.0), lambda(0.0), lambda_6(0.0), m0Squared(0.0), mHSquared(0.0), mH_in_GeV(0.0), resultFlag(0){};
	resultForOutput(const resultForOutput &o): cutoff_in_GeV(o.cutoff_in_GeV), vev(o.vev), y_t(o.y_t), y_b(o.y_b), lambda(o.lambda), lambda_6(o.lambda_6), 
	                 m0Squared(o.m0Squared), mHSquared(o.mHSquared), mH_in_GeV(o.mH_in_GeV), resultFlag(o.resultFlag){};
};

bool printResultsVectorToStream(const std::vector< resultForOutput > &results, std::ostream &output);

//counts the number of signchanges
int countNumberOfSignChanges(const std::map< double, double > &firstDerivative);






#endif
