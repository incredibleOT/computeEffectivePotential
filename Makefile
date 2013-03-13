CXX=g++

# CXXFLAGS=-DF_
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic 


OFILES = effectivePotential.o scanEffectivePotentialHelper.o


# ##########################################

scanEffectivePotentialByParameters: scanEffectivePotentialByParameters.o ${OFILES}
	${CXX} ${CXXFLAGS} -o $@ $^ ${LIBS}

scanEffectivePotentialByParameters.o: scanEffectivePotentialByParameters.cc ${OFILES}
	${CXX} ${CXXFLAGS} -c -o $@ $<


### ofiles
effectivePotential.o: effectivePotential.cc effectivePotential.h
	${CXX} ${CXXFLAGS} -c -o $@ $<

scanEffectivePotentialHelper.o: scanEffectivePotentialHelper.cc scanEffectivePotentialHelper.h
	${CXX} ${CXXFLAGS} -c -o $@ $<
	
#############################################

clean: 
	rm *.o *~