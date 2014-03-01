#ifndef _LOAD_FULL_PARTICLE_H_
#define _LOAD_FULL_PARTICLE_H_

struct full_particle {
  float pos[6], mass, energy;
  int32_t type;
  //Particle ID, Unique assigned internal halo id, 
  // Internal Halo ID, External Halo ID
  int64_t id, a_hid, hid, ehid;
};

void load_full_particles(char *filename, struct halo **h, int64_t *num_h, 
			 struct full_particle **p, int64_t *num_p, float *bnds);

#define GROUP_LIST h
#define parent n_core
#define RADIUS r
#define mm num_child_particles


#endif /* _LOAD_FULL_PARTICLE_H_ */
