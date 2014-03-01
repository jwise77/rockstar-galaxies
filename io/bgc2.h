/* From Cameron McBride, Copyright (c) 2011.
   Explicit permission granted to distribute with Rockstar under the GPLv3.
*/
#pragma once
#ifndef BGC_HEADER_DEFINED
#define BGC_HEADER_DEFINED 1

#include <inttypes.h>
#include <limits.h>

#if LONG_MAX < LLONG_MAX
#warning Compiling in 32-bit mode is not recommended; correct operation with datasets larger than 2GB is not guaranteed.
#endif

#define BGC_MAGIC ((uint64_t)0x1234567801010101ll)
#define BGC2_HEADER_SIZE 1024

/* for simplicity, we waste bits in the header and use a minimum of types 
 * in addition, we tried to stack all int types first followed by double */
typedef struct bgc2_header {
    uint64_t magic;             /* A magic number to identify this as a BGC file. */
    int64_t version;            /* File version number. */

    int64_t num_files;          /* number of files output is distributed into */
    int64_t file_id;            /* this files ID (number) if multiple files are output */
    int64_t snapshot;           /* Snapshot ID */

    int64_t format_group_data;  /* output group data format identifier, see enum gdata_format below */
    int64_t format_part_data;   /* output particle data format identifier, see enum pdata_format below */

    int64_t group_type;         /* FOF, SO, etc */
    int64_t ngroups;            /* number of groups stored LOCALLY in this file */
    int64_t ngroups_total;      /* number of groups stored GLOBALLY over all output BGC files */

    int64_t npart;              /* number of particles bound in groups LOCALLY in this file */
    int64_t npart_total;        /* number of particles bound in groups GLOBALLY over all output BGC files */
    int64_t npart_orig;         /* number of particles from original simulation input */

    int64_t max_npart;          /* maximum number of particles in one group LOCALLY in this file */
    int64_t max_npart_total;    /* maximum number of particles in one group GLOBALLY over all output BGC files */

    int64_t min_group_part;     /* minimum number of particles in a group */
    int64_t valid_part_ids;     /* valid particle IDs mean they match input snapshot */

    double linkinglength;       /* for FOF halos, what linking length is used */
    double overdensity;         /* mostly SO: overdensity with respect to mean */

    double time;                /* time of the input snapshot */
    double redshift;            /* redshift of the input snapshot */

    double box_size;            /* input BoxSize */
    double box_min[3];          /* alternative to center, Gadget assumes (0,0,0) */
    double bounds[6];           /* Spatial bounds of the halo centers contained in this file */

    double part_mass;           /* mass of particles if only one */

    double Omega0;              /* Matter density at z=0 */
    double OmegaLambda;         /* Dark energy density at z=0 */

    /* these define the units, but are NOT ALWAYS SET */
    double Hubble0;             /* in whatever units are used */
    double GravConst;           /* in whatever units are used */

    uint8_t padding[BGC2_HEADER_SIZE - ( 36 * 8 )];
} BGC2_HEADER;

enum gdata_format {
    GDATA_FORMAT_ID = 10,
    GDATA_FORMAT_RM = 20,
    GDATA_FORMAT_RMPV = 30,
    GDATA_FORMAT_RMPVMAX = 40,
};

enum group_type {
    GTYPE_FOF = 10,             /* standard 3D FOF */
    GTYPE_SOFOF = 20,           /* SO from an initial FOF halo */
    GTYPE_SO = 30,              /* SO using full DM, centering halos is an art. */
    GTYPE_SOB = 40,             /* SO "bound": using bound DM for mass */
};

typedef struct {
    int64_t id;
    int64_t parent_id;
    uint64_t npart;
} GROUP_DATA_ID;

typedef struct {
    int64_t id;
    int64_t parent_id;
    uint64_t npart;
    float radius;               /* halo radius. SO is well defined -- FOF not so much */
    float mass;
} GROUP_DATA_RM;

typedef struct {
    int64_t id;
    int64_t parent_id;
    uint64_t npart;
    float radius;
    float mass;
    float pos[3];
    float vel[3];
} GROUP_DATA_RMPV;

typedef struct {
    int64_t id;
    int64_t parent_id;
    uint64_t npart;
    uint64_t npart_self; /* Excluding substructure */
    float radius;
    float mass;
    float pos[3];
    float vel[3];
    float vmax, rvmax;
} GROUP_DATA_RMPVMAX;

/* Get the size of the GROUP DATA (gdata) structure based on FORMAT number */
#ifdef BGC2_SIZE
inline size_t
bgc_sizeof_gdata( const int64_t gdata_format )
{
    switch ( gdata_format ) {
        case GDATA_FORMAT_ID:
            return sizeof( GROUP_DATA_ID );
        case GDATA_FORMAT_RM:
            return sizeof( GROUP_DATA_RM );
        case GDATA_FORMAT_RMPV:
            return sizeof( GROUP_DATA_RMPV );
        case GDATA_FORMAT_RMPVMAX:
            return sizeof( GROUP_DATA_RMPVMAX );
    }

    fprintf( stderr, "ERROR: unknown group data format!  (format = %" PRId64 ")\n", gdata_format );
    return 0;
}
#endif /*BGC2_SIZE*/

enum pdata_format {
    PDATA_FORMAT_NULL = 0,
    PDATA_FORMAT_ID = 10,
    PDATA_FORMAT_IDBE = 15,
    PDATA_FORMAT_POS = 20,
    PDATA_FORMAT_POSBE = 25,
    PDATA_FORMAT_PV = 30,
    PDATA_FORMAT_PVBE = 35,
    PDATA_FORMAT_PVM = 40,
    PDATA_FORMAT_PVMBE = 45,
    PDATA_FORMAT_GPVM = 50
};

#ifdef BGC2_SIZE
inline int
bgc_format_includes_be( int64_t format_id )
{
    /* this MUST be kept in sync with the above enum defining the formats */
    return ( format_id % 10 ) == 5;
}
#endif /*BGC2_SIZE*/

typedef struct {
    int64_t part_id;
} PARTICLE_DATA_ID;

typedef struct {
    int64_t part_id;
    float binding_energy;
} PARTICLE_DATA_IDBE;

typedef struct {
    int64_t part_id;
    float pos[3];
} PARTICLE_DATA_POS;

typedef struct {
    int64_t part_id;
    float pos[3];
    float binding_energy;
} PARTICLE_DATA_POSBE;

typedef struct {
    int64_t part_id;
    float pos[3];
    float vel[3];
} PARTICLE_DATA_PV;

typedef struct {
    int64_t part_id;
    float pos[3];
    float vel[3];
    float binding_energy;
} PARTICLE_DATA_PVBE;

typedef struct {
    int64_t part_id;
    float pos[3];
    float vel[3];
    float mass;
} PARTICLE_DATA_PVM;

typedef struct {
    int64_t part_id;
    float pos[3];
    float vel[3];
    float mass;
    float binding_energy;
} PARTICLE_DATA_PVMBE;

typedef struct {
    int64_t part_id;
    int64_t group_id;
    float pos[3];
    float vel[3];
    float mass;
} PARTICLE_DATA_GPVM;

/* Get the size of the PARTICLE DATA (pdata) structure based on FORMAT number */
#ifdef BGC2_SIZE
inline size_t
bgc_sizeof_pdata( const int64_t pdata_format )
{
    switch ( pdata_format ) {
        case PDATA_FORMAT_NULL:
            return 0;
        case PDATA_FORMAT_ID:
            return sizeof( PARTICLE_DATA_ID );
        case PDATA_FORMAT_IDBE:
            return sizeof( PARTICLE_DATA_IDBE );
        case PDATA_FORMAT_POS:
            return sizeof( PARTICLE_DATA_POS );
        case PDATA_FORMAT_POSBE:
            return sizeof( PARTICLE_DATA_POSBE );
        case PDATA_FORMAT_PV:
            return sizeof( PARTICLE_DATA_PV );
        case PDATA_FORMAT_PVBE:
            return sizeof( PARTICLE_DATA_PVBE );
        case PDATA_FORMAT_PVM:
            return sizeof( PARTICLE_DATA_PVM );
        case PDATA_FORMAT_PVMBE:
            return sizeof( PARTICLE_DATA_PVMBE );
        case PDATA_FORMAT_GPVM:
            return sizeof( PARTICLE_DATA_GPVM );
    }

    fprintf( stderr, "ERROR: unknown particle data format!  (format = %" PRId64 ")\n",
             pdata_format );
    return 0;
}
#endif /*BGC2_SIZE*/

#endif

// vim: ts=4 sw=4 sts=4 expandtab
