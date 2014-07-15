ifdef HDF5
CPPFLAGS += -I${HDF5}/include -W -Wall -ansi
LDFLAGS += -L${HDF5}/lib -Wl,-rpath -Wl,${HDF5}/lib
endif

LIBS += -lhdf5_hl -lhdf5
CFLAGS += -O2 -g

all: fast5tofastq

fast5tofastq: fast5tofastq.c Makefile
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} $< -o $@ ${LIBS}

clean:
distclean:
	rm -f *~ fast5tofastq
