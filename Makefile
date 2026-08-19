CXX = g++
CXXFLAGS = -O3 -fopenmp
TARGET = hpmc
SRC = HP-MC.cpp

all: $(TARGET) setup

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

setup:
	bash setup.sh

clean:
	rm -f $(TARGET)
