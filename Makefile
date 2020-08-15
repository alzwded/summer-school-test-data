all: genplate.exe
	cd plate && rm -f *.txt *.csv && ../genplate 100 2000 100 10
	#cd plate && rm -f *.txt *.csv && ../genplate 16 2000 100 10

genplate.exe: genplate.cpp
	g++ -o genplate -O3 genplate.cpp -lm
