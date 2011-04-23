CC=g++
AR=ar
CFLAGS=-c -Wall -O2
LDFLAGS=

SOURCES=CoordGeodetic.cpp \
	CoordTopographic.cpp \
	Eci.cpp \
	Globals.cpp \
	Julian.cpp \
	Observer.cpp \
	SGP4.cpp \
	SatelliteOrbit.cpp \
        Timespan.cpp \
	Tle.cpp \
	Vector.cpp
OBJECTS=$(SOURCES:.cpp=.o)
SGP4LIB=libsgp4.a

TESTPROG=RunTest
TESTSOURCES=RunTest.cpp
TESTOBJECTS=$(TESTSOURCES:.cpp=.o)

all: $(SGP4LIB) ${TESTPROG}

${SGP4LIB}: ${OBJECTS}
	${AR} -rcs -o $@ ${OBJECTS}

${TESTPROG}: ${SGP4LIB} ${TESTOBJECTS}
	$(CC) ${TESTOBJECTS} $(LDFLAGS) -static -L. -lsgp4 -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o ${SGP4LIB} ${TESTPROG}
