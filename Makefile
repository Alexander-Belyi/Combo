# Various flags
CXX  = g++
LINK = $(CXX)
#CXXFLAGS = -I -Wall -g -O0 -DDEBUG -std=c++17
CXXFLAGS = -I -pedantic -Wall -Wextra -Wreorder -Wcast-align -Wcast-qual \
	-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self \
	-Wlogical-op -Wnoexcept \
	-Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow \
	-Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 \
	-Wswitch-default -Wundef -Wno-unused \
	-Wparentheses -Wconversion \
	-Werror \
	-O3 -funroll-loops -pipe -std=c++17
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




