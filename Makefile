
include Makefile.in

OBJECTS  = test.o C2Decomp.o Alloc.o TransposeX2Y.o TransposeY2Z.o TransposeZ2Y.o TransposeY2X.o MemSplitMerge.cpp IO.o

all: TEST

test.o: test.cpp C2Decomp.hpp 
	$(CC) $(CFLAGS) -c $< 

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

C2Decomp.o: C2Decomp.cpp C2Decomp.hpp
	$(CC) $(CFLAGS) -c $< 

TEST:  $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LIBF) 

clean: 
	rm -rf   *.o 


