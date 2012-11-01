// -*- Mode: c++ -*- 

#ifndef POINTPOT_H
#define POINTPOT_H

#include <stdio.h>
#include <stdint.h>

#define restrict __restrict__

namespace PointPot {

  class BinaryReader {
  private:
    FILE *fp;
    void xread(void *ptr, size_t size, size_t count);

  public:
    BinaryReader(FILE *fp);
    long read_integer();
    void check_magic(long);
    double read_double();
    void read_doubles(double *, size_t count);
    void read_integers(int32_t *, size_t count);
    void read_chars(char *, size_t count);
  };

  class BasePointPot {

  public:
    // setup locations and orientations of points
    void setup_current_points(const double pos_i[3], const double quat_i[4],
                              const double pos_j[3], const double quat_j[4]);

    void setup_current_points(const double pos_i[3], const double rot_i[3][3],
                              const double pos_j[3], const double rot_j[3][3]);

    void get_current_point_set_sizes(int &N_i, int &N_j);

    // applies virtual methods on current points and takes into account scaling
    double calculate_current_potential();
    void calculate_current_forces_torques(double *force_i, double *torque_i,
                                          double *force_j, double *torque_j);

    inline double get_r_max() {return 2*r_sphere + r_cutoff;}

  protected:
    // takes owernship of points, deleted in deconstructor
    BasePointPot(size_t N_base_points, const double *base_points_positions, const int *base_points_types,
                 int truncate_points, double r_sphere, double r_cutoff, double scale);
    virtual ~BasePointPot();

    virtual double calculate_a_potential(const int N_i, const double * restrict pos_i,
                                         const int * restrict types_i,
                                         const int N_j, const double * restrict pos_j,
                                         const int * restrict types_j) = 0;
    
    // force and torques on particle i
    virtual void calculate_a_force_torque(double *force, double *torque,
                                          const double com[3],
                                          const int N_i, const double * restrict pos_i,
                                          const int * restrict types_i,
                                          const int N_j, const double * restrict pos_j,
                                          const int * restrict types_j) = 0;

    inline double get_r_cutoff() {return r_cutoff;}

  private:
    const int N_base_points;
    const double * restrict base_point_positions;
    const int * restrict base_point_types;
    const int truncate_points;
    const double r_sphere;
    const double r_cutoff;
    const double scale;

    int N_cur_points_i, N_cur_points_j;
    double * restrict cur_point_positions_i, * restrict cur_point_positions_j;
    int * restrict cur_point_types_i, * restrict cur_point_types_j;
    double com_i[3], com_j[3];

    void setup_a_collection(const double pos[3], const double rot[3][3],
                            const double pos_other[3],
                            int &N, double * restrict &positions, int * restrict &types);
  };

  class CoulombPointPot : public BasePointPot {

  public:
    static CoulombPointPot *load(BinaryReader &br);

    double calculate_a_potential(const int N_i, const double *pos_i, const int *types_i,
                                 const int N_j, const double *pos_j, const int *types_j);
    
    void calculate_a_force_torque(double *force, double *torque, const double com[3],
                                  const int N_i, const double *pos_i, const int *types_i,
                                  const int N_j, const double *pos_j, const int *types_j);

  protected:
    CoulombPointPot(size_t N_base_points, const double *base_points_positions, const int *base_points_types,
                    int truncate_points, double r_sphere, double r_cutoff, double scale);


  };
}

#endif /* POINTPOT_H */
