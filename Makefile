all: compile run
compile: cudaLeapint.cu cudaLeapint.cuh
	nvcc -o cudaLeapint cudaLeapint.cu -lm -arch=sm_35 -rdc=true
run: compile
	./cudaLeapint 
clean: 
	rm cudaLeapint
	rm *.txt
	rm *.csv
	rm -rf SinglePoints/
