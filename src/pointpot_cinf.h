// -*- Mode: c++ -*- 

#ifndef POINTPOT_CINF_H
#define POINTPOT_CINF_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  void * POINTPOT_COUL_load(FILE *fp);
  void POINTPOT_COUL_free(void *);

  double POINTPOT_BASE_calculate_potential(void *ptr,
                                           const double pos_i[3], const double rot_i[3][3],
                                           const double pos_j[3], const double rot_j[3][3]);

  void POINTPOT_BASE_calculate_forces_torques(void *ptr,
                                              double *force_i, double *torque_i,
                                              double *force_j, double *torque_j,
                                              const double pos_i[3], const double rot_i[3][3],
                                              const double pos_j[3], const double rot_j[3][3]);

  void POINTPOT_BASE_calculate_point_charge_set_sizes(void *ptr, int *N_i, int *N_j,
                                                      const double pos_i[3], const double rot_i[3][3],
                                                      const double pos_j[3], const double rot_j[3][3]);


#ifdef __cplusplus
}
#endif

#endif /* POINTPOT_CINF_H */
