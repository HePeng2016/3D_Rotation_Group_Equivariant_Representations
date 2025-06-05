objects = include/rotation.o 
CC = g++  -std=c++11  -g  -I include 
all:  prepareall 
	$(CC)  -o rotation main.cpp $(objects) 
prepareall:  subsystem 
subsystem: 
	$(MAKE) -C  include/ 
clean:  cleansub 
	rm  rotation 
cleansub:
	rm $(objects)

