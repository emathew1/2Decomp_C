
include Makefile.in

OBJECTS  = test.o C2Decomp.o 

all: TEST

test.o: test.cpp C2Decomp.hpp 
	$(CC) $(CFLAGS) -c $< 

C2Decomp.o: C2Decomp.cpp C2Decomp.hpp 
	$(CC) $(CFLAGS) -c $< 

TEST:  $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LIBF) 

clean: 
	rm -rf   *.o 


