
include Makefile.in

TEST_OBJECTS  = test.o C2Decomp.o Alloc.o TransposeX2Y.o TransposeY2Z.o TransposeZ2Y.o TransposeY2X.o MemSplitMerge.cpp IO.o Best2DGrid.o Halo.o
TESTMAJORINDEXING_OBJECTS  = test_majorindexing.o C2Decomp.o Alloc.o TransposeX2Y.o TransposeY2Z.o TransposeZ2Y.o TransposeY2X.o MemSplitMerge.cpp IO.o Best2DGrid.o Halo.o
HALOTEST_OBJECTS  = halo_test.o C2Decomp.o Alloc.o TransposeX2Y.o TransposeY2Z.o TransposeZ2Y.o TransposeY2X.o MemSplitMerge.cpp IO.o Best2DGrid.o Halo.o

all: TEST TEST_MAJORINDEXING HALO_TEST

TransposeX2Y.o: TransposeX2Y.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

TransposeY2Z.o: TransposeY2Z.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

TransposeZ2Y.o: TransposeZ2Y.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

TransposeY2X.o: TransposeY2X.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

Alloc.o: Alloc.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

MemSplitMerge.o: MemSplitMerge.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

IO.o: IO.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

Best2DGrid.o: Best2DGrid.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

Halo.o: Halo.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

C2Decomp.o: C2Decomp.cpp C2Decomp.hpp
	$(CC) $(CFLAGS) -c $< 

test.o: test.cpp C2Decomp.hpp 
	$(CC) $(CFLAGS) -c $< 

halo_test.o: halo_test.cpp C2Decomp.hpp 
	$(CC) $(CFLAGS) -c $< 

TEST:  $(TEST_OBJECTS)
	$(CC) $(CFLAGS) $(TEST_OBJECTS) -o $@ $(LIBF) 

TEST_MAJORINDEXING:  $(TESTMAJORINDEXING_OBJECTS)
	$(CC) $(CFLAGS) $(TESTMAJORINDEXING_OBJECTS) -o $@ $(LIBF) 

HALO_TEST:  $(HALOTEST_OBJECTS)
	$(CC) $(CFLAGS) $(HALOTEST_OBJECTS) -o $@ $(LIBF) 



clean: 
	rm -rf   *.o 


