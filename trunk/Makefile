
CC=g++

#CCFLAGS=-g -Wall -fPIC
CCFLAGS  = -std=c++0x -g -O3 $(shell pkg-config --cflags cairo)

LDFLAGS = -Wl -shared  $(shell pkg-config --libs cairo)

all: color_compress0r.o 
	$(CC) -o color_compress0r.so $^ $(LDFLAGS)

%.o : %.cpp 
	$(CC) $(CCFLAGS) -c -o $@ $<

clean:
	- rm -f *.o *.so *~
