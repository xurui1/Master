LFLAGS = -lm  -O3 -fopenmp

main: main.cpp 
	g++ $(LFLAGS) -o $@ $(MOBLIB) main.cpp $(LIBS)
