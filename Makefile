all: compile distribution run plot
compile: cudaLeapint.cu cudaLeapint.cuh
	nvcc -o cudaLeapint cudaLeapint.cu -lm -arch=sm_35 -rdc=true
	gcc -o diskDist_diskInvader diskDist_diskInvader.c -lm -std=gnu99
distribution: diskDist_diskInvader
	./diskDist_diskInvader
run: compile initDiskGalaxy_forDarrensCudaNBody.csv
	./cudaLeapint initDiskGalaxy_forDarrensCudaNBody.csv 
clean: 
	rm cudaLeapint
	rm diskDist_diskInvader
	rm *.txt
	rm *.csv
	rm -rf SinglePoints/
plot: run initDiskGalaxy_forDarrensCudaNBodySim.txt
	./animateDensity.py initDiskGalaxy_forDarrensCudaNBodySim.txt
