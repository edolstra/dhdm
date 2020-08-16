dhdm: dhdm.cc
	g++ -std=c++17 -O3 -g dhdm.cc -o dhdm -pthread -lGLEW -lglfw -lGL -lpng -Wall -lfmt -losdCPU -lboost_iostreams
