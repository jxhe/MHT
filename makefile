main: cokus.o dataset.o mht-alpha.o mht-data.o mht-estimate.o mht-inference.o  mht-model.o strtokenizer.o utils.o
	g++ -O2 -fopenmp cokus.o dataset.o mht-alpha.o mht-data.o mht-estimate.o mht-inference.o  mht-model.o strtokenizer.o utils.o -o main

cokus.o:cokus.cpp cokus.h
	g++ -c -O2 -fopenmp cokus.cpp

dataset.o: dataset.cpp dataset.h
	g++ -c -O2 -fopenmp dataset.cpp

mht-alpha.o: mht-alpha.cpp mht-alpha.h
	g++ -c -O2 -fopenmp mht-alpha.cpp

mht-data.o: mht-data.cpp mht-data.h
	g++ -c -O2 -fopenmp mht-data.cpp

mht-estimate.o: mht-estimate.cpp mht-estimate.h
	g++ -c -O2 -fopenmp mht-estimate.cpp

mht-inference.o: mht-inference.cpp mht-inference.h
	g++ -c -O2 -fopenmp mht-inference.cpp

mht-model.o: mht-model.cpp mht-model.h
	g++ -c -O2 -fopenmp mht-model.cpp

strtokenizer.o: strtokenizer.cpp strtokenizer.h
	g++ -c -O2 -fopenmp strtokenizer.cpp

utils.o: utils.cpp utils.h
	g++ -c -O2 -fopenmp utils.cpp

clean:
	rm *.o