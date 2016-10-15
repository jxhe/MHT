main: cokus.o dataset.o lda-alpha.o lda-data.o lda-estimate.o lda-inference.o  lda-model.o strtokenizer.o utils.o
	g++ -O2 -fopenmp cokus.o dataset.o lda-alpha.o lda-data.o lda-estimate.o lda-inference.o  lda-model.o strtokenizer.o utils.o -o main

cokus.o:cokus.cpp cokus.h
	g++ -c -O2 -fopenmp cokus.cpp

dataset.o: dataset.cpp dataset.h
	g++ -c -O2 -fopenmp dataset.cpp

lda-alpha.o: lda-alpha.cpp lda-alpha.h
	g++ -c -O2 -fopenmp lda-alpha.cpp

lda-data.o: lda-data.cpp lda-data.h
	g++ -c -O2 -fopenmp lda-data.cpp

lda-estimate.o: lda-estimate.cpp lda-estimate.h
	g++ -c -O2 -fopenmp lda-estimate.cpp

lda-inference.o: lda-inference.cpp lda-inference.h
	g++ -c -O2 -fopenmp lda-inference.cpp

lda-model.o: lda-model.cpp lda-model.h
	g++ -c -O2 -fopenmp lda-model.cpp

strtokenizer.o: strtokenizer.cpp strtokenizer.h
	g++ -c -O2 -fopenmp strtokenizer.cpp

utils.o: utils.cpp utils.h
	g++ -c -O2 -fopenmp utils.cpp

clean:
	rm *.o