/* pointpot pair forcefield for lammps
 * Copyright (C) 2012 Matt Hagy <hagy@gatech.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "assert.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"

#include "pair_point_pot.h"


using namespace LAMMPS_NS;

#define dok3 for(int k=0; k<3; k++)

PairPointPot::PairPointPot(LAMMPS *lmp) : Pair(lmp)
{
  setflag = NULL;
  respa_enable = 0;
  point_pot = NULL;
}

PairPointPot::~PairPointPot()
{
  if (setflag) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
  }
  if (point_pot) delete point_pot;
}

void PairPointPot::allocate()
{
}

void PairPointPot::settings(int narg, char **arg)
{
  assert(point_pot == NULL);
  if (narg != 0) error->all("Illegal pair_style command for point_pot");
}

void PairPointPot::coeff(int narg, char **arg)
{
  int n = atom->ntypes;
  if (n != 1) error->all("Pair style currently only works with 1 atom type");
  if (narg != 3) error->all("Illegal pair_coeff command for point_pot");

  const char *filename = arg[2];

  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    char str[512];
    snprintf(str, sizeof(str), "Cannot open data file %.256s", filename);
    error->all(str);
  }
  PointPot::BinaryReader br = PointPot::BinaryReader(fp);
  point_pot = PointPot::CoulombPointPot::load(br);
  fclose(fp);

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 1;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  allocated = 1;
}

void PairPointPot::init_style()
{
  if (!atom->torque_flag || !(atom->quat_flag || atom->mu_flag))
    error->all("Pair point_force requires atom attributes torque and dipole or quaternion");
  int irequest = neighbor->request(this);
}

double PairPointPot::init_one(int i, int j)
{
  assert(point_pot != NULL);
  return point_pot->get_r_max();
}

void PairPointPot::write_restart(FILE *fp)
{
  error->all("cannot be written to restart");
}

void PairPointPot::write_restart_settings(FILE *fp)
{
  error->all("cannot be written to restart");
}

void PairPointPot::read_restart(FILE *fp)
{
  error->all("cannot be read from restart");
}

void PairPointPot::read_restart_settings(FILE *fp)
{
  error->all("cannot be read from restart");
}

static inline int
epsilon_eq(double a, double b, double epsilon)
{
  return fabs(a-b) <= epsilon;
}

static inline void
set_mat_identity(double mat[3][3])
{
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      mat[i][j] = i==j ? 1 : 0;
}

static inline void
canonical_rot_mat(const double *start, const double *end, double mat[3][3])
{
  const double theta_epsilon = 1e-8 * M_PI;
  double theta = acos(MathExtra::dot3(start, end));
  if (epsilon_eq(theta, 0, theta_epsilon)) {
    set_mat_identity(mat);
  } else if (epsilon_eq(theta, M_PI, theta_epsilon)) {
    set_mat_identity(mat);
    mat[2][2] = -1.0;
  } else {
    double axis[3];
    MathExtra::cross3(start, end, axis);
    double quat[4];
    MathExtra::axisangle_to_quat(axis, theta, quat);
    MathExtra::quat_to_mat(quat, mat);
  }
}

void PairPointPot::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  assert(point_pot != NULL);
  double r_max_sqr = point_pot->get_r_max() * point_pot->get_r_max();

  int n_evaled = 0;


  // loop over neighbors of my atoms
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];

    // could rotate points i into place ahead of time???
    double pos_i[3];
    dok3{ pos_i[k] = x[i][k]; }

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      double factor_lj;
      if (j < nall) factor_lj = 1.0;
      else {
        factor_lj = special_lj[j/nall];
        j %= nall;
      }

      double pos_j[3];
      dok3{ pos_j[k] = x[j][k]; }

      { // check if we're beyond cutoff distance
        double r2=0.0;
        dok3 {
          double x = pos_i[k] - pos_j[k];
          r2 += x*x;
        }
        //printf("r2=%.5g r_max_sqr=%.5g\n", r2, r_max_sqr);
        if (r2 >= r_max_sqr)
          continue;
      }

      n_evaled ++;

      if (atom->quat_flag) {
        double quat_i[4], quat_j[4];
        for (int k=0; k<4; k++) {
          quat_i[k] = atom->quat[i][k];
          quat_j[k] = atom->quat[j][k];
        }
        point_pot->setup_current_points(pos_i, quat_i, pos_j, quat_j);
      } else {
        double rot_i[3][3], rot_j[3][3];
        const double v_z[3] = {0.0, 0.0, 1.0};
        canonical_rot_mat(v_z, atom->mu[i], rot_i);
        canonical_rot_mat(v_z, atom->mu[j], rot_j);
        assert(atom->mu_flag);
        point_pot->setup_current_points(pos_i, rot_i, pos_j, rot_j);
      }

      double force_i[3], torque_i[3], force_j[3], torque_j[3];
      point_pot->calculate_current_forces_torques(force_i, torque_i,
                                                  force_j, torque_j);

      dok3 { 
        f[i][k] += factor_lj * force_i[k];
        tor[i][k] += factor_lj * torque_i[k];
      }

      if (newton_pair || j<nlocal) {
        dok3 { 
          f[j][k] += factor_lj * force_j[k];
          tor[j][k] += factor_lj * torque_j[k];
        }
      }

      if (evflag) {
        double r12[3];
        dok3 { r12[k] = pos_j[k] - pos_i[k]; }

        double evdwl = factor_lj * point_pot->calculate_current_potential();
      
        ev_tally_xyz(i,j,nlocal,newton_pair,
                     evdwl, 0.0,
                     factor_lj * force_i[0],
                     factor_lj * force_i[1],
                     factor_lj * force_i[2],
                     -r12[0], -r12[1], -r12[2]);
      }
    }
  }

  //printf("n_evaled=%d inum=%d r_max=%.4g \n", n_evaled, inum, point_pot->get_r_max());
  if (vflag_fdotr) virial_compute();
}
