#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include <assert.h>
#include "config_vars.h"
#include "universal_constants.h"
#include "potential.h"
#include "halo.h"

#define MAX_SCALE_BINS 50
#define MIN_PART_PER_BIN 15
#define MIN_SCALE_PART (100)

inline double c_to_f(double c) {
  double cp1 = 1.0+c;
  return (c*cp1 / (log1p(c)*cp1 - c));
}

double f_to_c(double f) {
  double c = f;
  double tc = c_to_f(c);
  double new_c;
  while (fabs((f-tc) / f) > 1e-7) {
    double tc2 = c_to_f(c+0.1);
    double slope = (tc2 - tc) / 0.1;
    new_c = c + (f-tc)/slope;
    if (new_c < 0) c/=2;
    else c = new_c;
    tc = c_to_f(c);
  }
  return c;
}

float estimate_scale_radius(float mvir, float rvir, float vmax, float rvmax, float scale)
{
  float f, c, vm2;
  if (!mvir || !rvir || !vmax) return (rvmax / RMAX_TO_RS);
  vm2 = vmax/VMAX_CONST;
  vm2 = vm2*vm2 * scale;
  f = (rvir/1.0e3)*vm2 / (mvir*RS_CONSTANT);
  if (f < 4.625) return (rvmax / RMAX_TO_RS);
  c = f_to_c(f);
  if (c <= 0) return (rvmax / RMAX_TO_RS);
  return (rvir/c);
}

static inline float nfw_menc(float r, float rs) {
  return (log((rs+r)/rs)-r/(rs+r));
}

float chi2_scale(float rs, float *bin_r, float *weights, int64_t num_bins) {
  int64_t i;
  float bin_mass = nfw_menc(bin_r[num_bins], rs)/((double)num_bins);
  float last_menc = 0;
  float chi2 = 0, dx;
  for (i=0; i<num_bins; i++) {
    float menc = nfw_menc(bin_r[i+1], rs);
    dx = weights[i]*((menc-last_menc) - bin_mass)/bin_mass;
    chi2 += dx*dx;
    last_menc = menc;
  }
  return chi2;
}

float calc_scale_from_bins(float rs, float *bin_r, float *weights, int64_t num_bins) {
  if (!(rs > bin_r[0]) || !(rs < bin_r[num_bins])) rs = bin_r[(num_bins+1)/2];
  float initial_rs = rs;
  float chi2 = chi2_scale(rs, bin_r, weights, num_bins);
  float drs = rs/10.0;
  int64_t iter = 0;
  float last_chi2 = 1e30;
  if (!(rs>0)) return initial_rs;
  while (fabs(chi2-last_chi2) > 0.005*chi2 && iter < MIN_PART_PER_BIN) {
    last_chi2 = chi2;
    if (rs+drs > bin_r[num_bins]) drs = 0.1*(bin_r[num_bins]-rs);
    if (rs-drs < bin_r[0]) drs = 0.1*(rs-bin_r[0]);
    float chi2_right = chi2_scale(rs+drs, bin_r, weights, num_bins);
    float chi2_left = chi2_scale(rs-drs, bin_r, weights, num_bins);
    float dx = 0.5*(chi2_right-chi2_left)/drs;
    float dx2 = (chi2_right+chi2_left-2.0*chi2)/(drs*drs);
    float move = 0;
    if (dx2 != 0) move = -dx/dx2;
    if (!move) return rs;
    if (rs+move > 4*rs) move = 3.0*rs;
    if (rs+move < 0.25*rs) move = -0.75*rs;
    float new_rs = rs+move;
    if (new_rs > bin_r[num_bins]) new_rs = rs + 0.8*(bin_r[num_bins]-rs);
    else if (new_rs < bin_r[0]) new_rs = rs + 0.8*(bin_r[0]-rs);
    float new_chi2 = chi2_scale(new_rs, bin_r, weights, num_bins);
    if (chi2_right < new_chi2) { new_chi2 = chi2_right; new_rs = rs+drs; }
    if (chi2_left < new_chi2) { new_chi2 = chi2_left; new_rs = rs-drs; }
    if (last_chi2 < new_chi2) { drs*=0.5; continue; }
    drs = fabs(rs-new_rs)/10.0;
    rs = new_rs;
    chi2 = new_chi2;
    iter++;
  }
  return rs;
}

void calc_scale_radius(struct halo *h, float mvir, float rvir, float vmax, float rvmax, float scale, struct potential *po, int64_t total_p, int64_t bound)
{
  float bin_r[MAX_SCALE_BINS+1]={0};
  float weights[MAX_SCALE_BINS]={0};
  int64_t i,num_bins=0,analyze_p=0,ppbin=0;
  float rs = h->klypin_rs = estimate_scale_radius(mvir, rvir, vmax, rvmax, scale);
  double total_mass = 0, mass, mpbin=0,max_mass=0;
  for (i=0; i<total_p-1; i++) {
     if (bound && (po[i].pe < po[i].ke)) continue;
     analyze_p++;
     total_mass += po[i].mass;
     if (po[i].mass > max_mass) max_mass = po[i].mass;
  }
  if (analyze_p < MIN_SCALE_PART) { h->rs = rs; return; }
  ppbin = ceil(((double)analyze_p/(double)MAX_SCALE_BINS));
  if (ppbin < MIN_PART_PER_BIN) 
    ppbin = (analyze_p) / ((int64_t)(analyze_p/MIN_PART_PER_BIN));
  num_bins = analyze_p / ppbin;
  if (num_bins>MAX_SCALE_BINS) num_bins=MAX_SCALE_BINS;
  mpbin = (double)total_mass / ((double)num_bins);
  if (mpbin < max_mass) mpbin = max_mass;
  num_bins=0;
  bin_r[0] = 0;
  mass = 0;
  for (i=0; i<total_p-1; i++) {
    if (bound && (po[i].pe < po[i].ke)) continue;
    mass += po[i].mass;
    if (mass >= mpbin) {
      num_bins++;
      float r1 = sqrt(po[i].r2), r2 = sqrt(po[i+1].r2);
      bin_r[num_bins] = 1e3*0.5*(r1+r2);
      mass -= mpbin;
      assert(po[i].mass > 0);
      float diff = mass/po[i].mass;
      if (i>0) { bin_r[num_bins]-=diff*(r1-sqrt(po[i-1].r2)); }
      weights[num_bins-1] = 1;
    }
  }
  
  for (i=0; i<num_bins; i++)
    if (bin_r[i]*1e-3 < 3*FORCE_RES) weights[i] = 0.1;

  h->rs = calc_scale_from_bins(rs, bin_r, weights, num_bins);
}
