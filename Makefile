all: plate cube

plate: genplate.exe
	cd plate && rm -f *.txt *.csv && ../genplate 100 2000 100 10
	#cd plate && rm -f *.txt *.csv && ../genplate 16 2000 100 10
.PHONY: plate

cube: gencube.exe
	#cd cube && rm -f *.txt *.csv && ../gencube 100 2000 100 10
	cd cube && rm -f *.txt *.csv && ../gencube 100 2000 100 16 
.PHONY: cube

%.exe: %.cpp
	g++ -fopenmp -o $@ -O3 $< -lm -lgomp
