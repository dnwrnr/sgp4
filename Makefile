CC=g++
AR=ar

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS=-c -Wall -O0 -g -pedantic -Wextra -Wconversion -march=native -Wno-long-long
else
	CFLAGS=-c -Wall -O2    -pedantic -Wextra -Wconversion -march=native -Wno-long-long
endif

LDFLAGS=

SOURCES=CoordGeodetic.cpp \
	CoordTopographic.cpp \
	Eci.cpp \
	Globals.cpp \
	Julian.cpp \
	Observer.cpp \
	OrbitalElements.cpp \
	SGP4.cpp \
	Timespan.cpp \
	Tle.cpp \
	Vector.cpp
OBJECTS=$(SOURCES:.cpp=.o)
SGP4LIB=libsgp4.a

TESTPROG=RunTest
TESTPROGSOURCES=RunTest.cpp
TESTPROGOBJECTS=$(TESTPROGSOURCES:.cpp=.o)

SATTRACK=SatTrack
SATTRACKSOURCES=SatTrack.cpp
SATTRACKOBJECTS=$(SATTRACKSOURCES:.cpp=.o)

PASSPREDICT=PassPredict
PASSPREDICTSOURCES=PassPredict.cpp
PASSPREDICTOBJECTS=$(PASSPREDICTSOURCES:.cpp=.o)

all: $(SGP4LIB) ${TESTPROG} ${SATTRACK} ${PASSPREDICT}

${SGP4LIB}: ${OBJECTS}
	${AR} -rcs $@ ${OBJECTS}

${TESTPROG}: ${SGP4LIB} ${TESTPROGOBJECTS}
	$(CC) ${TESTPROGOBJECTS} $(LDFLAGS) -L. -lsgp4 -o $@

${SATTRACK}: ${SGP4LIB} ${SATTRACKOBJECTS}
	${CC} ${SATTRACKOBJECTS} ${LDFLAGS} -L. -lsgp4 -o $@

${PASSPREDICT}: ${SGP4LIB} ${PASSPREDICTOBJECTS}
	${CC} ${PASSPREDICTOBJECTS} ${LDFLAGS} -L. -lsgp4 -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o ${SGP4LIB} ${TESTPROG} ${SATTRACK} ${PASSPREDICT}
