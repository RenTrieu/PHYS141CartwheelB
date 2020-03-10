all: compile distribution run plot
compile: cudaLeapint.cu cudaLeapint.cuh
	nvcc -o cudaLeapint cudaLeapint.cu -lm -arch=sm_35 -rdc=true
	gcc -o diskDist_diskInvader diskDist_diskInvader.c -lm -std=gnu99
distribution: diskDist_diskInvader
	./diskDist_diskInvader
massPoint: 
	cat initDiskGalaxy_forDarrensCudaNBody.csv | sed '10002,13312d' > temp.csv
	cat temp.csv | sed '9216, 10000d' > massInvader.csv
run: compile distribution massPoint 
	./cudaLeapint initDiskGalaxy_forDarrensCudaNBody.csv 
	./cudaLeapint massInvader.csv
clean: 
	rm cudaLeapint
	rm diskDist_diskInvader
	rm *.txt
	rm *.csv
	rm -rf SinglePoints/
plot: run initDiskGalaxy_forDarrensCudaNBodySim.txt
	./animateDensity.py initDiskGalaxy_forDarrensCudaNBodySim.txt
	./animateDensity.py massInvaderSim.txt

