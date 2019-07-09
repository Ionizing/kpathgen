#pragma once

#include "poscar.h"
#include "matops.h"
#include "error.h"

#include "spglib/spglib.h"
#include "spglib/spg_database.h"

/*
 * Formatting related stuff
 * From POSCAR struct to Cellp and versa
 *
 *  Cellp  ->  Cell plus
 *  Add comment line, type_names members
 */
typedef struct {
  int     size;
  int     n_types;
  double (*lattice)[3]; /* 3x3 matrix col major */
  int     *types;
  double (*position)[3];
  char    comment[1 << 10];
  char   (*type_names)[10];
} Cellp; 

Cellp* cellp_cell_alloc(const int n_atoms, const int n_types);
ERROR_CODE cellp_cell_free(Cellp* cell);
ERROR_CODE cellp_cell_set(Cellp* cell,
                          const char   comment[],
                          const double lattice[3][3],
                          const int    types[],
                          const double position[][3],
                          const int    n_atoms,
                          const char   type_names[][10],
                          const int    n_types);
ERROR_CODE cellp_cell_print(const Cellp* cell);


Cellp* psc_poscar2spg_cell(const POSCAR* poscar);
POSCAR* spg_cell2psc_poscar(const Cellp* cell);




/*
 *  
 */

// number of atoms is returned
int kpt_standardize_cell(   Cellp*   cell,
                        const bool   to_primitive, 
                        const bool   no_ideal, 
                        const double prec);

// space group number is returned, [1, 230]
int kpt_get_point_group(const Cellp* cell,
                        const double prec,
                              char   symbol[7]);

// hall number is returned, [1, 530]
int kpt_get_hall_number(const Cellp* cell,
                        const double prec);
