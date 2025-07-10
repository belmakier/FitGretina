INSTALLDIR=$(HOME)/.local/

CC = g++
CFLAGS = `root-config --cflags` -O3 -g -fPIC
ROOTLIBS = `root-config --libs --glibs` -lMathMore -lMinuit2

all : libFitGretina SingleGretinaFit pParGretinaFit

libFitGretina : src/FitGretina.cc src/FitGretina.hh
	$(CC) $(CFLAGS) -shared -o lib/libFitGretina.so src/FitGretina.cc $(ROOTLIBS) $(LIBS)

SingleGretinaFit : src/SingleFit.cc
	$(CC) $(CFLAGS) -o bin/SingleGretinaFit src/SingleFit.cc $(ROOTLIBS) $(LIBS) -L./lib/ -lFitGretina

pParGretinaFit : src/pParFit.cc
	$(CC) $(CFLAGS) -o bin/pParGretinaFit src/pParFit.cc $(ROOTLIBS) $(LIBS) -L./lib/ -lFitGretina

install :
	cp lib/libFitGretina.so ${INSTALLDIR}/lib/
	cp bin/* ${INSTALLDIR}/bin/

clean :
	rm lib/libFitGretina.so 
	rm bin/* 
