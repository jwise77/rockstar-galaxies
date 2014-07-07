#ifndef HALO_H
#define HALO_H

#define HALO_FORMAT_REVISION 1

#include <stdint.h>

struct halo {
  int64_t id;
  float pos[6], corevel[3], bulkvel[3];
  float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms,
    J[3], energy, spin, alt_m[4], Xoff, Voff, b_to_a, c_to_a, A[3],
    b_to_a2, c_to_a2, A2[3],
    bullock_spin, kin_to_pot, m_pe_b, m_pe_d;
  int64_t num_p, num_child_particles, p_start, desc, flags, n_core;
  float min_pos_err, min_vel_err, min_bulkvel_err;

  int32_t type;
  float sm, gas, bh, peak_density, av_density;
  //float half_sm_radius, sm_vrms 
};

struct extra_halo_info {
  int64_t child, next_cochild, prev_cochild;
  int64_t sub_of, ph;
  float max_metric, volume;
  double x_orth_matrix[3][3], v_orth_matrix[3][3];
  double x_eig[3], v_eig[3];
};

#endif /* HALO_H */
