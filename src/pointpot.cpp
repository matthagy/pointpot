
#include <assert.h>
#include <math.h>

#include "debug.h"
#include "pointpot.h"
#include "math_extra.h"

using namespace PointPot;


#define dok3 for(int k=0; k<3; k++)


/* * * * * * *
 * Binary IO *
 * * * * * * */

BinaryReader::BinaryReader(FILE *fp) :
  fp(fp)
{
  assert(fp != NULL);
}

void BinaryReader::xread(void *ptr, size_t size, size_t count)
{
  size_t r = fread(ptr, size, count, fp);
  if (r != count) {
    Fatal("failed read %lu for size %lu count %lu", r, size, count);
  }
}

long BinaryReader::read_integer()
{
  int32_t i;
  xread(&i, sizeof(int32_t), 1);
  return (long)i;
}

void BinaryReader::check_magic(long magic)
{
  long i = read_integer();
  if (i != magic) {
    Fatal("bad magic %lX; expected %lX", i, magic);
  }
}  
 
double BinaryReader::read_double()
{
  double d;
  xread(&d, sizeof(double), 1);
  return d;
}

void BinaryReader::read_integers(int32_t *ptr, size_t count)
{
  xread(ptr, sizeof(int32_t), count);
}

void BinaryReader::read_doubles(double *ptr, size_t count)
{
  xread(ptr, sizeof(double), count);
}

void BinaryReader::read_chars(char *ptr, size_t count)
{
  xread(ptr, sizeof(char), count);
}

/* * * * * * * * * * * * * * *
 * Consturctor & Destructor  *
 * * * * * * * * * * * * * * */

BasePointPot::BasePointPot(size_t N_base_points, const double *base_point_positions, 
                                const int *base_point_types,
                                int truncate_points, double r_sphere, double r_cutoff, double scale) :
  N_base_points(N_base_points), base_point_positions(base_point_positions),
  base_point_types(base_point_types), truncate_points(truncate_points),
  r_sphere(r_sphere), r_cutoff(r_cutoff), scale(scale)
{
  assert(N_base_points > 0);
  assert(base_point_positions != 0);
  assert(base_point_types != 0);
  assert(r_sphere > 0);
  assert(r_cutoff > 0);
  assert(scale != 0);

  cur_point_positions_i = new double[3*N_base_points];
  cur_point_positions_j = new double[3*N_base_points];

  if (!truncate_points) {
    N_cur_points_i = N_cur_points_j = N_base_points;
    cur_point_types_i = cur_point_types_j = (int *)base_point_types;
  } else {
    N_cur_points_i = N_cur_points_j = 0;
    cur_point_types_i = new int[N_base_points];
    cur_point_types_j = new int[N_base_points];
  }
}

BasePointPot::~BasePointPot()
{
  delete []base_point_positions;
  delete []base_point_types;
  delete []cur_point_positions_i;
  delete []cur_point_positions_j;
  if (truncate_points) {
    delete [] cur_point_types_i;
    delete [] cur_point_types_j;
  }
}


/* * * * * * * * * * * * * * * *
 * Setup of Point Collections  *
 * * * * * * * * * * * * * * * */

void BasePointPot::setup_current_points(const double pos_i[3], const double quat_i[4],
                                        const double pos_j[3], const double quat_j[4])
{
  double rot_i[3][3], rot_j[3][3];
  MathExtra::quat_to_mat(quat_i, rot_i);
  MathExtra::quat_to_mat(quat_j, rot_j);
  setup_current_points(pos_i, rot_i, pos_j, rot_j);
}

void BasePointPot::setup_current_points(const double pos_i[3], const double rot_i[3][3],
                                        const double pos_j[3], const double rot_j[3][3])
{
  dok3{com_i[k] = pos_i[k];}
  dok3{com_j[k] = pos_j[k];}
  setup_a_collection(pos_i, rot_i, pos_j, N_cur_points_i, cur_point_positions_i, cur_point_types_i);
  setup_a_collection(pos_j, rot_j, pos_i, N_cur_points_j, cur_point_positions_j, cur_point_types_j);
}

void BasePointPot::get_current_point_set_sizes(int &N_i, int &N_j)
{
  N_i = N_cur_points_i;
  N_j = N_cur_points_j;
}

static inline void
show_vec(const char *name, double scale, const double vec[3])
{
  printf("%.64s=[%.4f, %.4f, %.4f]\n", name, 
         scale*vec[0], scale*vec[1], scale*vec[2]);
}

static inline void
calculate_dividing_plane(double plane_normal[3], double &plane_d,
                         const double pos[3], const double pos_other[3], double r_sphere)
{
  double r[3];
  dok3{ r[k] = pos[k] - pos_other[k]; }
  MathExtra::normalize3(r, plane_normal);
    
  double nearest_other_point[3];
  double scale = r_sphere / sqrt(MathExtra::dot3(r, r));
  dok3{ nearest_other_point[k] = pos_other[k] + r[k] * scale; }

  plane_d = -MathExtra::dot3(plane_normal, nearest_other_point);
}


void BasePointPot::setup_a_collection(const double pos[3], const double rot[3][3],
                                      const double pos_other[3],
                                      int &N, double * restrict &positions, int * restrict &types)
{
  if (!truncate_points) {
    assert(N == N_base_points);
    assert(types == base_point_types);
    for (int i=0; i<N; i++) {
      int offset = 3*i;
      double rotated[3];
      MathExtra::times_column3(rot, base_point_positions + offset, rotated);
      dok3{ positions[offset + k] = rotated[k] + pos[k]; }
    }
  } else {
    N = 0;
    assert(types != base_point_types);

    double plane_normal[3], plane_d;

    calculate_dividing_plane(plane_normal, plane_d, pos, pos_other, r_sphere);

    //printf("plane [%.3f %.3f %.3f] %.3g\n", plane_normal[0], plane_normal[1], plane_normal[2], plane_d);

    for (int i=0; i<N_base_points; i++) {

      double transformed[3];
      MathExtra::times_column3(rot, base_point_positions + 3*i, transformed);
      dok3{ transformed[k] += pos[k]; }

      double plane_distance = MathExtra::dot3(plane_normal, transformed) + plane_d;
      assert(plane_distance > 0);
      if (plane_distance <= r_cutoff) {
        dok3{ positions[N*3 + k] = transformed[k]; }
        types[N] = base_point_types[i];
        N += 1;
      }
    }
    //printf("setup %d of %d\n", N, N_base_points);
    assert(N >= 0);
    assert(N <= N_base_points);
  }
}

/* * * * * * * * * * * * * * * * * *
 * Calculations on current setups  *
 * * * * * * * * * * * * * * * * * */

double BasePointPot::calculate_current_potential()
{
  return scale * calculate_a_potential(N_cur_points_i, cur_point_positions_i, cur_point_types_i,
                                       N_cur_points_j, cur_point_positions_j, cur_point_types_j);
}

void BasePointPot::calculate_current_forces_torques(double *force_i, double *torque_i,
                                                    double *force_j, double *torque_j)
{
  //printf("N_cur_points_i=%d N_cur_points_j=%d\n", N_cur_points_i, N_cur_points_j);
  calculate_a_force_torque(force_i, torque_i, com_i,
                           N_cur_points_i, cur_point_positions_i, cur_point_types_i,
                           N_cur_points_j, cur_point_positions_j, cur_point_types_j);
  calculate_a_force_torque(force_j, torque_j, com_j,
                           N_cur_points_j, cur_point_positions_j, cur_point_types_j,
                           N_cur_points_i, cur_point_positions_i, cur_point_types_i);
  dok3 {
    force_i[k] *= scale;
    torque_i[k] *= scale;
    force_j[k] *= scale;
    torque_j[k] *= scale;
  }
}


/* * * * * * * * * * * * * *
 * Coulomb Point Potential *
 * * * * * * * * * * * * * */

double CoulombPointPot::calculate_a_potential(const int N_i, const double * restrict pos_i, const int * restrict types_i,
                                              const int N_j, const double * restrict pos_j, const int * restrict types_j)
{
  double r_cutoff_sqr = get_r_cutoff() * get_r_cutoff();
  double acc_U = 0.0;
  int n_shift = 0;

  // nasty, but efficient loops
  const double *pos_i_ptr=pos_i, *pos_i_end=pos_i + 3*N_i;
  const int *type_i_ptr = types_i;
  while (pos_i_ptr != pos_i_end) {
    double x = *(pos_i_ptr++);
    double y = *(pos_i_ptr++);
    double z = *(pos_i_ptr++);
    int type_i = *(type_i_ptr++);

    const double *pos_j_ptr=pos_j, *pos_j_end=pos_j + 3*N_j;
    const int *type_j_ptr = types_j;
    while (pos_j_ptr != pos_j_end) {
      double dx = x - *(pos_j_ptr++);
      double dy = y - *(pos_j_ptr++);
      double dz = z - *(pos_j_ptr++);
      int same_type = type_i == *(type_j_ptr++);

      double r_sqr = dx*dx + dy*dy + dz*dz;
      if (r_sqr < r_cutoff_sqr) {
        int s = (same_type<<1) - 1; //1 if same_type else -1
        n_shift += s;
        acc_U += (double)s / sqrt(r_sqr);
      }
    }
  }

  acc_U -= n_shift * (1.0 / get_r_cutoff());
  return acc_U;
}

void CoulombPointPot::calculate_a_force_torque(double *com_force, double *com_torque,
                                               const double com[3],
                                               const int N_i, const double * restrict pos_i, 
                                               const int * restrict types_i,
                                               const int N_j, const double * restrict pos_j,
                                               const int * restrict types_j)
{
  //printf("N_i=%d N_j=%d\n", N_i, N_j);

  // zero accumulators
  dok3{ com_force[k] = com_torque[k] = 0.0; }

  // iterate over pairs
  double r_cutoff_sqr = get_r_cutoff() * get_r_cutoff();

  for (int i=0; i<N_i; i++) {
    double loc[3];
    dok3{ loc[k] = pos_i[3*i + k]; }
    int tpi = types_i[i];
    double acc_force[3] = {0.0, 0.0, 0.0};
                
    for (int j=0; j<N_j; j++) {
      double r[3];
      dok3{ r[k] = loc[k] - pos_j[3*j + k];}
      double r2 = MathExtra::dot3(r, r);

      if (r2 < r_cutoff_sqr) {
        double f_l = (tpi == types_j[j] ? 1.0 : -1.0) / (r2 * sqrt(r2));
        dok3{ acc_force[k] += r[k] * f_l; }
      }
    }

    // accumulate translational force
    dok3{ com_force[k] += acc_force[k]; }

    // torque arm
    double r[3];
    dok3{ r[k] = loc[k] - com[k]; }

    double tau[3];
    MathExtra::cross3(r, acc_force, tau);

    // accumulate torque
    dok3{ com_torque[k] += tau[k]; }
  }

}

CoulombPointPot::CoulombPointPot(size_t N_base_points, const double *base_points_positions, const int *base_points_types,
                                 int truncate_points, double r_sphere, double r_cutoff, double scale) :
  BasePointPot(N_base_points, base_points_positions, base_points_types,
               truncate_points, r_sphere, r_cutoff, scale)
{
}


CoulombPointPot * CoulombPointPot::load(BinaryReader &br)
{
  const long magic = 0x734234;
  br.check_magic(magic);
  int N_base_points = br.read_integer();
  double *base_point_positions = new double[3*N_base_points];
  br.read_doubles(base_point_positions, 3*N_base_points);

  int32_t *tmp_types = new int32_t[N_base_points];
  br.read_integers(tmp_types, N_base_points);  
  int *base_point_types = new int[N_base_points];
  for (int i=0; i<N_base_points; i++) {
    int value = (int)tmp_types[i];
    assert(value == 1 || value == -1);
    base_point_types[i] = value;
  }
  delete tmp_types;

  int truncate_points = br.read_integer();
  double r_sphere = br.read_double();
  double r_cutoff = br.read_double();
  double scale = br.read_double();

  br.check_magic(magic);
  
  return new CoulombPointPot(N_base_points, base_point_positions, base_point_types,
                             truncate_points, r_sphere, r_cutoff, scale);
}
