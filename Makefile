CC=g++
CFLAGS=-c -Wall -g
LDFLAGS=
SOURCES=Coord.cpp \
	Globals.cpp \
	Julian.cpp \
	main.cpp \
	SGDP4.cpp \
	Tle.cpp \
	Vector.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=sat

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)
