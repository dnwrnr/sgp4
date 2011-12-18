CC=g++
AR=ar

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS=-Wall -O0 -g -pedantic -Werror -Wextra -Wconversion
else
	CFLAGS=-Wall -O2    -pedantic -Werror -Wextra -Wconversion -DNDEBUG
endif

LDFLAGS=

SOURCES=Eci.cpp \
	Globals.cpp \
	Julian.cpp \
	Observer.cpp \
	OrbitalElements.cpp \
	SGP4.cpp \
	SolarPosition.cpp \
	Timespan.cpp \
	Tle.cpp
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
	${AR} rcs $@ ${OBJECTS}

${TESTPROG}: ${SGP4LIB} ${TESTPROGOBJECTS}
	$(CC) ${LDFLAGS} ${CFLAGS} ${TESTPROGOBJECTS} -static -L. -lsgp4 -o $@

${SATTRACK}: ${SGP4LIB} ${SATTRACKOBJECTS}
	${CC} ${LDFLAGS} ${CFLAGS} ${SATTRACKOBJECTS} -static -L. -lsgp4 -o $@

${PASSPREDICT}: ${SGP4LIB} ${PASSPREDICTOBJECTS}
	${CC} ${LDFLAGS} ${CFLAGS} ${PASSPREDICTOBJECTS} -static -L. -lsgp4 -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -rf *.o ${SGP4LIB} ${TESTPROG} ${SATTRACK} ${PASSPREDICT}
