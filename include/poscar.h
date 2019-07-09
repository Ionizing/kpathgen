#pragma once
#include "base.h"
#include "error.h"
#include "strops.h"
#include "matops.h"

/*
 * [[a_x, a_y, a_z],
 *  [b_x, b_y, b_z],
 *  [c_x, c_y, c_z]]
 */

typedef struct {
  char    comment[1 << 10];
  double  lattice[3][3];      // lattice vectors, row major
  int     n_types;            // number of types
  int    (*types);            // atom number of each type
  char   (*type_names)[10];   // atom type name
  int     n_atoms;            // total atoms
  double (*positions)[3];     // atom positions, stored in FRACTIONAL COORDINATES
} POSCAR;

POSCAR*      psc_poscar_alloc(const int n_atoms,
                              const int n_types);
ERROR_CODE   psc_poscar_free(POSCAR* poscar);
ERROR_CODE   psc_poscar_set(POSCAR* poscar,
                            const char*  comment,
                            const double lattice[3][3],
                            const int types[],
                            const char type_names[][10],
                            const int n_types,
                            const double positions[][3],
                            const int n_atoms);
POSCAR*      psc_poscar_read(const char* filename);
ERROR_CODE   psc_poscar_print(const POSCAR* poscar);
ERROR_CODE   psc_poscar_write(const POSCAR* poscar,
                              const char* filename);


