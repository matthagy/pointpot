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

#ifdef PAIR_CLASS

PairStyle(point_pot, PairPointPot)

#else

#ifndef LMP_PAIR_POINT_FORCE_H
#define LMP_PAIR_POINT_FORCE_H

#include "pair.h"
#include "stdio.h"

#include <pointpot.h>

namespace LAMMPS_NS {

class PairPointPot : public Pair {
 public:
  PairPointPot(LAMMPS *lmp);
  virtual ~PairPointPot();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  void allocate();
  PointPot::CoulombPointPot *point_pot;
};

 
}
#endif
#endif
