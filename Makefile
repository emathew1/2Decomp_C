
include Makefile.in

OBJECTS  = test.o C2Decomp.o Alloc.o TransposeX2Y.o MemSplitMerge.cpp

all: TEST

test.o: test.cpp C2Decomp.hpp 
	$(CC) $(CFLAGS) -c $< 

TransposeX2Y.o: TransposeX2Y.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

Alloc.o: Alloc.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

MemSplitMerge.o: MemSplitMerge.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

C2Decomp.o: C2Decomp.cpp C2Decomp.hpp
	$(CC) $(CFLAGS) -c $< 

TEST:  $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LIBF) 

clean: 
	rm -rf   *.o 


