#ifndef _NFW_H_
#define _NFW_H_

void calc_scale_radius(struct halo *h, float mvir, float rvir, float vmax, float rvmax, float scale, struct potential *po, int64_t total_p, int64_t bound);

#endif /* _NFW_H_ */
