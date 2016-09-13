#ifndef _IO_TIPSY_H_
#define _IO_TIPSY_H_
#define TIPSY_MAXDIM 3

#include "../particle.h"
#include <rpc/types.h>
#include <rpc/xdr.h>

struct tipsy_gas_particle {
    float mass;
    float pos[TIPSY_MAXDIM];
    float vel[TIPSY_MAXDIM];
    float rho;
    float temp;
    float hsmooth;
    float metals ;
    float phi ;
} ;

struct tipsy_dark_particle {
    float mass;
    float pos[TIPSY_MAXDIM];
    float vel[TIPSY_MAXDIM];
    float eps;
    float phi ;
} ;

struct tipsy_star_particle {
    float mass;
    float pos[TIPSY_MAXDIM];
    float vel[TIPSY_MAXDIM];
    float metals ;
    float tform ;
    float eps;
    float phi ;
} ;

struct tipsy_dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;

struct tipsy_dump_pad {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
    int pad;
} ;

struct tipsy_nc_dump {
  int    magic;
  double time;
  int    iHighWord;
  int    nbodies;
  int    ndim;
  int    code;

} ;

void load_particles_tipsy(char *filename, struct particle **p, int64_t *num_p);
int load_ids_tipsy(char *filename, struct tipsy_dump header, int **iords);
int tipsy_xdr_header(XDR *pxdrs,struct tipsy_dump *ph);
int tipsy_xdr_gas(XDR *pxdrs,struct tipsy_gas_particle *ph);
int tipsy_xdr_dark(XDR *pxdrs,struct tipsy_dark_particle *ph);
int tipsy_xdr_star(XDR *pxdrs,struct tipsy_star_particle *ph);

#endif /* _IO_TIPSY_H */

