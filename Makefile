all: compile run
compile: cudaLeapint.cu cudaLeapint.cuh
	nvcc -o cudaLeapint cudaLeapint.cu -lm -arch=sm_35 -rdc=true
	g++ -o nBodyLeapint nBodyLeapint.cpp OctTree.cpp
run: compile
	./cudaLeapint 
clean: 
	rm cudaLeapint

