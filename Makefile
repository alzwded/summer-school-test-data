all: plate cube16 cube32 cube

plate: genplate.exe
	mkdir -p plate
	cd plate && rm -f *.txt *.csv *.tec && ../genplate.exe 100 2000 100 100
.PHONY: plate

cube: gencube.exe
	mkdir -p cube
	cd cube && rm -f *.txt *.csv *.tec && ../gencube.exe 100 2000 100 16
.PHONY: cube

cube32: gencube.exe
	mkdir -p cube32
	cd cube32 && rm -f *.txt *.csv *.tec && ../gencube.exe 32 2000 100 16
.PHONY: cube32

cube16: gencube.exe
	mkdir -p cube16
	cd cube16 && rm -f *.txt *.csv *.tec && ../gencube.exe 16 2000 100 16 
.PHONY: cube16

%.exe: %.cpp
	g++ -march=sandybridge -fopenmp -o $@ -O3 $< -lm

clean:
	rm -rf *.exe cube cube16 cube32 plate *.obj *.ilk *.pdb
