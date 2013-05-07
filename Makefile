CXX=g++

# CXXFLAGS=-DF_
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic 


OFILES= effectivePotential.o scanEffectivePotentialHelper.o

effPotFiles=  effectivePotential_computePropagatorSum.cc effectivePotential_getFermionicContribution.cc effectivePotential.h

INCL=-I/opt/products/gsl/1.15/include

# ##########################################

scanEffectivePotentialByParameters: scanEffectivePotentialByParameters.o ${OFILES}
	${CXX} ${CXXFLAGS} -o $@ $^ ${LIBS}

scanEffectivePotentialByParameters.o: scanEffectivePotentialByParameters.cc ${OFILES}
	${CXX} ${CXXFLAGS} -c -o $@ $<

determineBoundAndPlotPotential: determineBoundAndPlotPotential.o effectivePotential.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LIBS}

determineBoundAndPlotPotential.o: determineBoundAndPlotPotential.cc effectivePotential.o
	${CXX} ${CXXFLAGS} -c -o $@ $<



### ofiles
effectivePotential.o: effectivePotential.cc   ${effPotFiles}
	${CXX} ${CXXFLAGS} -c -o $@ $< ${effPotOfiles}


scanEffectivePotentialHelper.o: scanEffectivePotentialHelper.cc scanEffectivePotentialHelper.h
	${CXX} ${CXXFLAGS} -c -o $@ $<
	
#############################################

clean: 
	rm *.o *~