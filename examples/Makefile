CFLAGS=-m64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -lm -fno-math-errno -Wall
OFLAGS = -O3 -std=c99
DEBUGFLAGS = -g -O3 -std=c99
PROFFLAGS = -g -pg -O2 -std=c99
CC = gcc

all:
	@make reg EXTRA_FLAGS="$(OFLAGS)"

debug:
	@make reg EXTRA_FLAGS="$(DEBUGFLAGS)"

prof:
	@make reg EXTRA_FLAGS="$(PROFFLAGS)"

reg:
	$(CC) -DCALC_POTENTIALS $(CFLAGS) bound_particle_assignments.c load_full_particles.c ../check_syscalls.c  ../io/stringparse.c ../io/io_util.c ../io/io_nchilada.c ../hubble.c ../config_vars.c ../potential.c -o bound_particle_assignments  $(EXTRA_FLAGS)
	$(CC) -DCALC_POTENTIALS $(CFLAGS) gen_grp_stats.c load_full_particles.c ../check_syscalls.c  ../io/stringparse.c ../io/io_util.c ../io/io_nchilada.c ../hubble.c ../config_vars.c ../potential.c -o gen_grp_stats  $(EXTRA_FLAGS)
	$(CC) $(CFLAGS) calc_bgc2_shapes.c load_bgc2.c ../check_syscalls.c ../io/io_util.c ../distance.c ../rockstar.c ../config_vars.c ../jacobi.c ../fun_times.c ../io/meta_io.c ../io/io_bgc2.c ../io/io_ascii.c ../io/stringparse.c ../io/io_art.c ../io/io_gadget.c ../io/io_tipsy.c ../potential.c ../bounds.c ../client.c ../config.c ../fof.c ../hubble.c ../integrate.c ../interleaving.c ../inthash.c ../merger.c ../nfw.c ../server.c ../subhalo_metric.c ../universe_time.c ../inet/address.c ../inet/rsocket.c ../inet/socket.c ../io/io_generic.c ../io/io_internal.c ../io/io_nchilada.c ../io/read_config.c -o calc_bgc2_shapes  $(EXTRA_FLAGS)
	$(CC) $(CFLAGS) bgc2_to_ascii_particles.c load_bgc2.c ../check_syscalls.c ../io/io_util.c -o bgc2_to_ascii_particles  $(EXTRA_FLAGS)
	$(CC) -DTEST_LOADFP $(CFLAGS) load_full_particles.c ../check_syscalls.c   ../io/stringparse.c ../config_vars.c -o load_full_particles  $(EXTRA_FLAGS)
	$(CC) -DCALC_POTENTIALS $(CFLAGS) calc_potentials.c load_full_particles.c ../check_syscalls.c  ../hubble.c ../io/stringparse.c ../config_vars.c ../potential.c -o calc_potentials  $(EXTRA_FLAGS)

clean:
	rm -f *~ 

