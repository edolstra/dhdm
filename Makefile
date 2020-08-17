dhdm: main.o dhdm.o mesh.o dsf.o subdivide.o diff.o render.o
	g++ -o $@ $^ -std=c++17 -O3 -g -pthread -lGLEW -lglfw -lGL -lpng -Wall -lfmt -losdCPU -lboost_iostreams

%.o: %.cc dhdm.hh mesh.hh diff.hh
	g++ -o $@ -c $< -std=c++17 -O3 -g

clean:
	rm -f *.o dhdm
