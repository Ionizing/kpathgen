#include "../include/spg_wrap.h"
#include "../include/error.h"

#define __DEBUG_ON
#include "../include/debug.h"


Cellp* cellp_cell_alloc(const int n_atoms, const int n_types) {
  Cellp* cell = NULL;
  if (n_atoms < 1 || n_types < 1) {
    error_puts("Bad Cell_plus alloc request: n_atoms or n_types < 1.");
    return NULL;
  } else {  }

  cell = (Cellp*)malloc(sizeof(Cellp));
  if (NULL == cell) {
    error_puts("Cannot allocate memory for Cellp object.");
    return NULL;
  } else {  }

  cell->size     = n_atoms;
  cell->n_types  = n_types;

  cell->lattice = (double (*)[3])malloc(sizeof(double[3][3]));
  if (NULL == cell->lattice) {
    error_puts("Cannot allocate memory for Cellp object.");

    free(cell);
    cell = NULL;

    return NULL;
  } else {  }

  cell->types = (int*) calloc(n_atoms, sizeof(int));
  if (NULL == cell->types) {
    error_puts("Cannot allocate memory for Cellp object");

    free(cell->lattice);
    cell->lattice = NULL;
    
    free(cell);
    cell = NULL;

    return NULL;
  } else {  }

  cell->position = (double (*)[3])calloc(n_atoms, sizeof(double[3]));
  if (NULL == cell->position) {
    error_puts("Cannot allocate memory for Cellp object.");

    free(cell->types);
    cell->types = NULL;

    free(cell->lattice);
    cell->lattice = NULL;

    free(cell);
    cell = NULL;

    return NULL;
  } else {  }

  cell->type_names = (char (*)[10])calloc(n_types, sizeof(char[10]));
  if (NULL == cell->type_names) {
    error_puts("Cannot allocate memory for Cellp object.");

    free(cell->position);
    cell->position = NULL;

    free(cell->types);
    cell->types = NULL;

    free(cell->lattice);
    cell->lattice = NULL;

    free(cell);
    cell = NULL;

    return NULL;
  } else {  }

  cell->comment[0] = '\0';
  return cell;
} // end of cellp_cell_alloc



ERROR_CODE cellp_cell_free(Cellp* cell) {
  if (NULL != cell) {
    if (NULL != cell->lattice) {
      free(cell->lattice);
      cell->lattice = NULL;
    } else {  }

    if (NULL != cell->position) {
      free(cell->position);
      cell->position = NULL;
    } else {  }

    if (NULL != cell->types) {
      free(cell->types);
      cell->types = NULL;
    }

    if (NULL != cell->type_names) {
      free(cell->type_names);
      cell->type_names = NULL;
    } else {  }

    free(cell);
  } else {  }

  return SUCCESS;
} // end of cellp_cell_free



ERROR_CODE cellp_cell_set(Cellp* cell,
                          const char   comment[],
                          const double lattice[3][3],
                          const int    types[],
                          const double position[][3],
                          const int    n_atoms,
                          const char   type_names[][10],
                          const int    n_types) {
  if (NULL == cell) {
    error_puts("NULL Cellp object pointer passed in when setting it.");
    return CELLP_NULL_OBJECT;
  } else {  }

  if (n_atoms != cell->size || n_types != cell->n_types) {
    error_puts("Bad Cellp object setting: n_atoms or n_types not consistent with provided cell.");
    return CELLP_INVALID_NATOMS_OR_NTYPES;
  } else {  }

  strncpy(cell->comment,   comment,    sizeof(cell->comment));
  memcpy(cell->lattice,    lattice,    sizeof(double[3][3]));
  memcpy(cell->types,      types,      sizeof(int) * n_atoms);
  memcpy(cell->position,   position,   sizeof(double[3]) * n_atoms);
  memcpy(cell->type_names, type_names, sizeof(char[10]) * n_types);

  mat_mat33_trans(cell->lattice);

  return SUCCESS;
}


Cellp* psc_poscar2spg_cell(const POSCAR* poscar) {
  if (NULL == poscar) {
    error_puts("NULL poscar object passed in when converting POSCAR to Cell.");
    return NULL;
  } else {  }

  const int n_atoms = poscar->n_atoms;
  const int n_types = poscar->n_types;

  Cellp* cell = cellp_cell_alloc(n_atoms, n_types);
  if (NULL == cell) {
    error_puts("Cannot allocate memory for cell object.");
    return NULL;
  } else {  }

 /* Converting [3, 4] to [0, 0, 0, 1, 1, 1, 1] */
  {
    int cnt = 0;
    for (int i=0; i!=poscar->n_types; ++i) {
      for (int j=0; j!=poscar->types[i]; ++j) {
        cell->types[cnt] = i;
        ++cnt;
      }
    }
  }

  cellp_cell_set(cell, poscar->comment, poscar->lattice, cell->types, poscar->positions,
      poscar->n_atoms, poscar->type_names, poscar->n_types);

  return cell;
}


POSCAR* spg_cell2psc_poscar(const Cellp* cell) {
  if (NULL == cell) {
    error_puts("NULL cell object passed in when converting Cell to POSCAR.");
    return NULL;
  } else {  }

  const int n_atoms = cell->size;
  const int n_types = cell->n_types;
  POSCAR* poscar = psc_poscar_alloc(n_atoms, n_types);
  if (NULL == poscar) {
    error_puts("Cannot allocate memory for poscar object.");
    return NULL;
  } else {  }

  int* types = (int*) calloc(n_types, sizeof(int));
  memset(types, sizeof(int) * n_types, 0);

  {
    int cnt = 0;
    for (int i=0; i!=n_atoms; ++i) {
      if (0 == i) {
        ++cnt;
      } else if (cell->types[i] != cell->types[i - 1]) {
        ++cnt;
      } else {  }
      ++types[cnt - 1];
    }
  }

  psc_poscar_set(poscar, 
      cell->comment, 
      cell->lattice,
      types, 
      cell->type_names, 
      n_types, 
      cell->position, 
      n_atoms);

  mat_mat33_trans(poscar->lattice);

  free(types);
  return poscar;
}


int kpt_standardize_cell(Cellp* cell,
                        const bool to_primitive,
                        const bool no_ideal,
                        const double prec) {
  if (NULL == cell) {
    error_puts("NULL cell object passed in when standardize cell.");
    return 0;
  } else {  }

  int res = spg_standardize_cell(cell->lattice, cell->position, cell->types, 
      cell->size,to_primitive, no_ideal, prec);
  return res;
}



int kpt_get_point_group(const Cellp* cell,
                        const double prec,
                        char symbol[7]) {
  return spg_get_schoenflies(symbol,
                             cell->lattice,
                             cell->position,
                             cell->types,
                             cell->size,
                             prec);
}


int kpt_get_hall_number(const Cellp* cell,
                        const double prec) {
  const int max_size = 50;
  int rotation[max_size][3][3];
  double translation[max_size][3];
  const int num_operations = spg_get_symmetry(rotation, 
                                              translation, 
                                              max_size,
                                              cell->lattice, 
                                              cell->position, 
                                              cell->types, 
                                              cell->size, 
                                              prec);

  return spg_get_hall_number_from_symmetry(rotation, 
                                           translation, 
                                           num_operations, 
                                           prec);
}
