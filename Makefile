# Various flags
CXX  = g++
LINK = $(CXX)
#CXXFLAGS = -I -Wall -g -O0 -DDEBUG -std=c++17
CXXFLAGS = -I -Wall -Wextra -Wreorder -O3 -funroll-loops -pipe -std=c++17
LFLAGS = -lm

TARGET  = comboCPP

HEADER  = Graph.h Combo.h
FILES = Graph.cc Combo.cc Main.cc

OBJECTS = $(FILES:.cc=.o)

$(TARGET): ${OBJECTS}
	$(LINK) $^ $(LFLAGS) -o $@

all: $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean:
	rm -f $(OBJECTS) $(TARGET)

# Compile and dependency
$(OBJECTS): $(HEADER) Makefile




