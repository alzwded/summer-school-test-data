all: plate cube16 cube32 cube

plate: genplate.exe
	mkdir -p plate
	cd plate && rm -f *.txt *.csv *.tec && ../genplate 100 2000 100 100
.PHONY: plate

cube: gencube.exe
	mkdir -p cube
	cd cube && rm -f *.txt *.csv *.tec && ../gencube 100 2000 100 16
.PHONY: cube

cube32: gencube.exe
	mkdir -p cube32
	cd cube32 && rm -f *.txt *.csv *.tec && ../gencube 32 2000 100 16
.PHONY: cube32

cube16: gencube.exe
	cd cube16 && rm -f *.txt *.csv *.tec && ../gencube 16 2000 100 16 
.PHONY: cube16

%.exe: %.cpp
	g++ -fopenmp -o $@ -O3 $< -lm

clean:
	rm -f *.exe cube cube16 cube32 plate *.obj
