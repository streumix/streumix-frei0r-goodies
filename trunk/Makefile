
# where to find frei0r.hpp ?
FREI0R_DIR = ../../../frei0r
# the plugins
PLUGINS= color_compress0r color_distanc0r selective_c0nv

CC=g++

CXXFLAGS  = -std=c++0x -g -O3 -I$(FREI0R_DIR)/include $(shell pkg-config --cflags cairo)
LDFLAGS = -lfrei0r_common $(shell pkg-config --libs cairo) -Wl,-shared 

all: $(PLUGINS)

frei0r_common.a: %.a : %.o
	ar rc $@ $<

$(PLUGINS) : % : %.cpp
	$(CC) $(CXXFLAGS) -c -o $@.so $< $(LDFLAGS)

clean:
	- rm -f *.o *.so *~
	
.PHONY: clean all

### depends
frei0r_common.o : frei0r_common.cpp frei0r_common.hpp 
$(PLUGINS) : frei0r_common.a frei0r_common.hpp hist_plot.hpp 
