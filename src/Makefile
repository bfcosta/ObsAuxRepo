CXX=em++
CXXFLAGS=-I/usr/include/eigen3/ -I../include -O3 -s WASM=1 -s NO_EXIT_RUNTIME=1 -s ALLOW_MEMORY_GROWTH=1 -s SINGLE_FILE=1
LDFLAGS=../lib/animation.o ../lib/bvhreader.o ../lib/kf2dinterpolator.o ../lib/kfinterpolator.o	../lib/kfselector.o ../lib/quad.o ../lib/rbf.o ../lib/skanalizer.o ../lib/skeleton.o
EXPFLAGS=-s EXPORTED_FUNCTIONS='["_getBvhProjection","_getNumberOfFrames","_getRbfSmoothedProjection","_getImprovedProjection","_getPoseError","_getSkfSkeleton","_getMocapSkeleton"]' -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]'
CC=emcc
CFLAGS=-I../include

all:	../js/projection.js ../js/multiproj.js

../js/projection.js: projection.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(EXPFLAGS) $(LDFLAGS)

../js/multiproj.js: multiproj.cpp ../lib/phelper.o
	$(CXX) -o $@ $< $(CXXFLAGS) $(EXPFLAGS) $(LDFLAGS) ../lib/phelper.o

../lib/quad.o: quad.c ../include/gluos.h ../include/gluint.h
	$(CC) -o $@ -c $< $(CFLAGS)

../lib/phelper.o: phelper.cpp ../include/phelper.h
	$(CC) -o $@ -c $< $(CFLAGS)
