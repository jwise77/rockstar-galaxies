#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "io_art.h"
#include "io_util.h"
#include "io_gadget.h"
#include "../universal_constants.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"

#define ART_BUFFER 100000

void art_detect_endianness(FILE *input, char *filename) {
  int32_t first_word, swapped;
  SWAP_ENDIANNESS=0;
  check_fread(&first_word, sizeof(int32_t), 1, input);
  check_limited_funread(&first_word, sizeof(int32_t), 1);
  if (first_word != ART_TEXT_SIZE) {
    SWAP_ENDIANNESS = 1;
    swapped = first_word;
    swap_endian_4byte((int8_t *)(&first_word));
  }
  if (first_word != ART_TEXT_SIZE) {
    fprintf(stderr, "[Error] Unrecognized ART file type in %s!\n", filename);
    fprintf(stderr, "Expected title size of %d; got %d (or %d if byte-swapped).\n", ART_TEXT_SIZE, first_word, swapped);
    exit(1);
  }
}

void art_detect_variant(FILE *input, char *filename) {
  int32_t first_word;
  check_fread(&first_word, sizeof(int32_t), 1, input);
  check_limited_funread(&first_word, sizeof(int32_t), 1);
  if (SWAP_ENDIANNESS) swap_endian_4byte((int8_t *)(&first_word));
  if (first_word == sizeof(struct art_header1)) return;
  if (first_word == sizeof(struct art_header1a)) { ART_VARIANT=1; return; }
  fprintf(stderr, "[Error] Unrecognized ART format in %s!\n", filename);
  fprintf(stderr, "[Error] Expected first header to be %d or %d bytes; got %d.\n", 
	  (int32_t)sizeof(struct art_header1), (int32_t)sizeof(struct art_header1a), first_word);
  exit(1);
}

void art_process_particles(struct particle *p, int64_t num_p, struct art_particle *ap, int64_t num_read, double NGRIDC)
{
  int64_t i, j;
  double Xscale = BOX_SIZE/NGRIDC;
  double Vscale = BOX_SIZE*100.0/(NGRIDC*SCALE_NOW);

  for (i=num_p; i<num_p+num_read; i++) {
    if (ART_VARIANT==0) {
      p[i].id = ap[i-num_p].id;
      if (SWAP_ENDIANNESS) swap_4byte_to_8byte((int32_t *)(void *)(&(p[i].id)));
      memcpy(p[i].pos, ap[i-num_p].pos, sizeof(float)*6);
    }
    for (j=0; j<3; j++) {
      p[i].pos[j]-=1.0;
      p[i].pos[j]*=Xscale;
      p[i].pos[j+3]*=Vscale;
    }
  }
}

void art_extract_header_info(struct art_header1 *p1)
{
  SCALE_NOW = p1->AEXPN;
  BOX_SIZE = p1->Box;
  Om = p1->Om0;
  Ol = p1->Oml0;
  h0 = p1->hubble;
  AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));
}

void _load_particles_art_v2_a(FILE *input, struct art_header1 *old_header)
{
  struct art_header1a part_art_header = {0};
  struct art_header2a part_art_header2 = {0};
  struct art_header3 part_art_header3 = {0};

  fread_fortran(&part_art_header, sizeof(struct art_header1a), 1, input, SWAP_ENDIANNESS);
  memcpy(old_header, &part_art_header, sizeof(struct art_header1));
  fread_fortran(&part_art_header2, sizeof(struct art_header2a), 1, input, SWAP_ENDIANNESS);
  fread_fortran(&part_art_header3, sizeof(struct art_header3), 1, input, SWAP_ENDIANNESS);
  PARTICLE_MASS = part_art_header.MassOne;
  TRIM_OVERLAP = part_art_header2.dBuffer;
  ROUND_AFTER_TRIM = 1.0/16.0;
}

void _load_particles_art_v2_b(FILE *input, struct particle **p, int64_t *num_p, uint32_t np, char *filename) {
  uint32_t p_in_record, total_read=0;
  char *buffer = NULL;
  int64_t i, bsize = 0;

  while ((total_read < np) && !feof(input)) {
    fread_fortran(&p_in_record, sizeof(uint32_t), 1, input, SWAP_ENDIANNESS);
    assert(total_read + p_in_record <= np);
    assert(p_in_record > 0);
    check_realloc_var(buffer, sizeof(float), bsize, 6*p_in_record);

    //Positions, velocities, then skip particle masses before reading IDs
    fread_fortran(buffer, sizeof(float), 3*p_in_record, input, SWAP_ENDIANNESS);
    fread_fortran(buffer + sizeof(float)*3*p_in_record, 
		  sizeof(float), 3*p_in_record, input, SWAP_ENDIANNESS);
    float *fbuffer = (float *)buffer;
    struct particle *tp = p[0]+(*num_p)+total_read;
    for (int64_t j=0; j<6; j++,fbuffer+=p_in_record)
      for (i=0; i<p_in_record; i++) tp[i].pos[j] = fbuffer[i];

    fread_fortran(buffer, sizeof(float), 3*p_in_record, input, SWAP_ENDIANNESS);
    int64_t *ibuffer = (int64_t *)(&(buffer[sizeof(float)*p_in_record]));
    for (i=0; i<p_in_record; i++) {
      if (SWAP_ENDIANNESS) swap_4byte_to_8byte((int32_t *)(&(ibuffer[i])));
      tp[i].id = ibuffer[i];
    }
    total_read += p_in_record;
  }
  *num_p = (*num_p) + np;
  fclose(input);
  if (buffer) free(buffer);
}


void load_particles_art(char *filename, struct particle **p, int64_t *num_p)
{
  FILE *input;
  uint32_t np, p_in_record, dummy;
  uint32_t total_read = 0, num_read;
  char title[ART_TEXT_SIZE];
  struct art_particle *pbuffer = NULL;
  struct art_header1 part_art_header = {0};
  struct art_header2 part_art_header2 = {0};
  struct art_header3 part_art_header3 = {0};

  input = check_fopen(filename, "r");
  art_detect_endianness(input, filename);
  fread_fortran(title, 1, ART_TEXT_SIZE, input, SWAP_ENDIANNESS);
  art_detect_variant(input, filename);

  if (ART_VARIANT==1) { _load_particles_art_v2_a(input, &part_art_header); }
  else {
    fread_fortran(&part_art_header, sizeof(struct art_header1), 1, input, SWAP_ENDIANNESS);
    fread_fortran(&part_art_header2, sizeof(struct art_header2), 1, input, SWAP_ENDIANNESS);
    fread_fortran(&part_art_header3, sizeof(struct art_header3), 1, input, SWAP_ENDIANNESS);
  }
  art_extract_header_info(&part_art_header);
  if (fabs(1.0-(Om+Ol)) > 1e-4) {
    fprintf(stderr, "[Error] Either cosmology is not flat or ART header structure has changed.\n");
    exit(1);
  }

  fread_fortran(&np, sizeof(uint32_t), 1, input, SWAP_ENDIANNESS);
  *p = check_realloc(*p, sizeof(struct particle)*((*num_p) + np),
		     "Allocating particles.");

  if (ART_VARIANT == 1) { _load_particles_art_v2_b(input, p, num_p, np, filename); return; }

  pbuffer = check_realloc(pbuffer, sizeof(struct art_particle)*ART_BUFFER,
			  "Allocating particle buffer.");

  while ((total_read < np) && !feof(input)) {
    fread_fortran(&p_in_record, sizeof(uint32_t), 1, input, SWAP_ENDIANNESS);
    check_fread(&dummy, sizeof(uint32_t), 1, input);
    while ((p_in_record > 0) && !feof(input)) {
      num_read = (p_in_record > ART_BUFFER) ? ART_BUFFER : p_in_record;
      if (SWAP_ENDIANNESS)
	num_read = fread_swap(pbuffer, sizeof(struct art_particle),
			      num_read, input);
      else num_read = check_fread(pbuffer, sizeof(struct art_particle),
			    num_read, input);
      art_process_particles(*p,*num_p,pbuffer,num_read,part_art_header.NGRIDC);
      p_in_record -= num_read;
      total_read += num_read;
      *num_p += num_read;
    }
    check_fread(&dummy, sizeof(uint32_t), 1, input);
  }
  free(pbuffer);
  fclose(input);
}
