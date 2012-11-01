
#include "pointpot.h"
#include "pointpot_cinf.h"

using namespace PointPot;

void * POINTPOT_COUL_load(FILE *fp)
{
  BinaryReader br = BinaryReader(fp);
  return (void *) CoulombPointPot::load(br);
}

void POINTPOT_COUL_free(void *ptr)
{
  delete reinterpret_cast<CoulombPointPot *>(ptr);
}

double POINTPOT_BASE_calculate_potential(void *ptr,
                                         const double pos_i[3], const double rot_i[3][3],
                                         const double pos_j[3], const double rot_j[3][3])
{
  BasePointPot *cpp = static_cast<BasePointPot *>(ptr);
  cpp->setup_current_points(pos_i, rot_i, pos_j, rot_j);
  return cpp->calculate_current_potential();
}


void POINTPOT_BASE_calculate_forces_torques(void *ptr,
                                    double *force_i, double *torque_i,
                                    double *force_j, double *torque_j,
                                    const double pos_i[3], const double rot_i[3][3],
                                    const double pos_j[3], const double rot_j[3][3])
{
  BasePointPot *cpp = static_cast<BasePointPot *>(ptr);
  cpp->setup_current_points(pos_i, rot_i, pos_j, rot_j);
  cpp->calculate_current_forces_torques(force_i, torque_i,
                                        force_j, torque_j);
}

void POINTPOT_BASE_calculate_point_charge_set_sizes(void *ptr, int *N_i, int *N_j,
                                                    const double pos_i[3], const double rot_i[3][3],
                                                    const double pos_j[3], const double rot_j[3][3])
{
  BasePointPot *cpp = static_cast<BasePointPot *>(ptr);
  cpp->setup_current_points(pos_i, rot_i, pos_j, rot_j);
  cpp->get_current_point_set_sizes(*N_i, *N_j);
}
