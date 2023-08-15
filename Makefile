CFLAGS=-m64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_DEFAULT_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_DARWIN_C_SOURCE -Wall -fno-math-errno -fPIC -I/usr/include/tirpc
LDFLAGS=-shared
OFLAGS = -lm -ltirpc -O3 -std=c99
DEBUGFLAGS = -lm -g -O3 -std=c99 #-Dinline= 
PROFFLAGS = -lm -g -pg -O2 -std=c99
CC = gcc
#CFILES = rockstar.c check_syscalls.c fof.c groupies.c subhalo_metric.c potential.c nfw.c jacobi.c fun_times.c interleaving.c universe_time.c hubble.c integrate.c distance.c config_vars.c config.c bounds.c inthash.c io/read_config.c client.c server.c merger.c inet/socket.c inet/rsocket.c inet/address.c io/meta_io.c io/io_internal.c io/io_ascii.c io/stringparse.c io/io_gadget.c io/io_generic.c io/io_art.c io/io_nchilada.c io/io_tipsy.c io/io_bgc2.c io/io_util.c io/io_arepo.c io/io_hdf5.c io/io_enzo.c
OFILES = main.o rockstar.o check_syscalls.o fof.o groupies.o subhalo_metric.o potential.o nfw.o jacobi.o fun_times.o interleaving.o universe_time.o hubble.o integrate.o distance.o config_vars.o config.o bounds.o inthash.o client.o server.o merger.o io/read_config.o inet/socket.o inet/rsocket.o inet/address.o io/meta_io.o io/io_internal.o io/io_ascii.o io/stringparse.o io/io_gadget.o io/io_generic.o io/io_art.o io/io_nchilada.o io/io_tipsy.o io/io_bgc2.o io/io_util.o io/io_arepo.o io/io_hdf5.o io/io_enzo.o
DIST_FLAGS =
HDF5_FLAGS = -DH5_USE_16_API -lhdf5 -DENABLE_HDF5 -I${HDF5_LIB}/include -L${HDF5_LIB}/lib
HDF5 = 1

ifeq ($(HDF5),1)
OFLAGS += $(HDF5_FLAGS)
endif

all: rockstar-galaxies

rockstar-galaxies: $(OFILES)
	$(CC) $(CFLAGS) $(OFILES) -o rockstar-galaxies $(OFLAGS) $(EXTRA_FLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $(OFLAGS) -c $*.c -o $@

lib: all
	$(CC) $(CFLAGS) $(LDFLAGS) $(OFILES) -o librockstar-galaxies.so $(OFLAGS) $(EXTRA_FLAGS)

#all:
#	@make reg EXTRA_FLAGS="$(OFLAGS)"

# with_hdf5:
# 	@make reg EXTRA_FLAGS="$(OFLAGS) $(HDF5_FLAGS)"

# debug:
# 	@make reg EXTRA_FLAGS="$(DEBUGFLAGS)"

# with_hdf5_debug:
# 	@make reg EXTRA_FLAGS="$(DEBUGFLAGS) $(HDF5_FLAGS)"

# prof:
# 	@make reg EXTRA_FLAGS="$(PROFFLAGS)"

# .REMAKE:

# dist: .REMAKE
# 	cd ../ ; perl -ne 'print "$$1\n" if (/VERSION\s*\"([^\"]+)/)' Rockstar-Galaxies/version.h > Rockstar-Galaxies/VERSION; tar -czvf rockstar-galaxies.tar.gz Rockstar-Galaxies/Makefile Rockstar-Galaxies/*.[ch] Rockstar-Galaxies/examples/Makefile Rockstar-Galaxies/[^sto]*/*.[ch] Rockstar-Galaxies/quickstart.cfg Rockstar-Galaxies/parallel.cfg Rockstar-Galaxies/scripts/*.pbs Rockstar-Galaxies/scripts/*.cfg Rockstar-Galaxies/scripts/*.pl Rockstar-Galaxies/SOURCE_LAYOUT Rockstar-Galaxies/README.pdf Rockstar-Galaxies/README Rockstar-Galaxies/LICENSE Rockstar-Galaxies/VERSION Rockstar-Galaxies/ACKNOWLEDGMENTS Rockstar-Galaxies/CHANGELOG; mv rockstar-galaxies.tar.gz Rockstar-Galaxies

# versiondist:
# 	$(MAKE) dist DIST_FLAGS="$(DIST_FLAGS)"
# 	rm -rf dist
# 	mkdir dist
# 	cd dist; tar xzf ../rockstar-galaxies.tar.gz ; perl -ne '/\#define.*VERSION\D*([\d\.rcRC-]+)/ && print $$1' Rockstar-Galaxies/version.h > NUMBER ; mv Rockstar-Galaxies Rockstar-Galaxies-`cat NUMBER`; tar czf rockstar-galaxies-`cat NUMBER`.tar.gz Rockstar-Galaxies-`cat NUMBER`

# reg:
# 	$(CC) $(CFLAGS) main.c $(CFILES) -o rockstar-galaxies  $(EXTRA_FLAGS)

# lib:
# 	$(CC) $(CFLAGS) $(LDFLAGS) $(CFILES) -o librockstar-galaxies.so $(OFLAGS)

# lib_with_hdf5:
# 	$(CC) $(CFLAGS) $(LDFLAGS) $(CFILES) -o librockstar-galaxies.so $(OFLAGS)  $(HDF5_FLAGS)

# bgc2:
# 	$(CC) $(CFLAGS) io/extra_bgc2.c util/redo_bgc2.c $(CFILES) -o util/finish_bgc2  $(OFLAGS)
# 	$(CC) $(CFLAGS) io/extra_bgc2.c util/bgc2_to_ascii.c $(CFILES) -o util/bgc2_to_ascii  $(OFLAGS)

# parents:
# 	$(CC) $(CFLAGS) util/find_parents.c io/stringparse.c check_syscalls.c  -o util/find_parents $(OFLAGS)

# substats:
# 	$(CC) $(CFLAGS) util/subhalo_stats.c $(CFILES) -o util/subhalo_stats  $(OFLAGS)


clean:
	rm -f *~ io/*~ inet/*~ util/*~ rockstar-galaxies util/redo_bgc2 util/subhalo_stats *.o io/*.o

