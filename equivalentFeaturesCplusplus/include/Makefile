objects = rotation.o
CC = g++  -g
MATCHFLAGS =  -I  wigner 
install : $(objects)
	
rotation.o  : rotation.cpp  rotation.h    
	$(CC)   -c  $<    $(MATCHFLAGS)
 
clean :
	rm  $(objects)
