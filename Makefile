LFLAGS = -lm -std=c++11 -O2 -fopenmp

main: main.cpp 
	g++ $(LFLAGS) -o $@ $(MOBLIB) main.cpp $(LIBS)
