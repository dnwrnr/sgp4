CC=g++
CFLAGS=-c -Wall -g
LDFLAGS=
SOURCES=Coord.cpp \
	Eci.cpp \
	Globals.cpp \
	Julian.cpp \
	Observer.cpp \
	SGDP4.cpp \
	Tle.cpp \
	Vector.cpp  \
	main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=sat

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)
