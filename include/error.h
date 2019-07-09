#pragma once

typedef enum {
  SUCCESS = 0,

  FILE_OPEN_FAILED,
  FILE_NULL_POINTER,

  POSCAR_ALLOC_FAIL,
  POSCAR_INVALID_ATOM_NUMBER,
  POSCAR_OPEN_FILE_FAIL,
  POSCAR_INVALID_SCALE,
  POSCAR_INCOMPLETE_CONTENT,
  POSCAR_NO_ATOM_TYPE_TAG,
  POSCAR_NO_ATOM_QUANTITY_LINE,
  POSCAR_ELEMENT_TAG_TOO_LONG,
  POSCAR_TYPES_NUMBERS_INCONSISTENT,
  POSCAR_NO_COORDINATE_TYPE_TAG,
  POSCAR_NULL_OBJECT,

  CELLP_NULL_OBJECT,
  CELLP_INVALID_NATOMS_OR_NTYPES,

  STRING_INVALID_VEC_SIZE,
  STRING_INVALID_STR_TO_FLOAT,
  STRING_INVALID_STR_TO_INT,
  STRING_NULL_OR_EMPTY,

  MATRIX_SINGULAR_MATRIX,
  MATRIX_INVALID_NROW,
} ERROR_CODE;

// Error handling
#define error_printf(fmt,...) fprintf(stderr, __BOLDRED "*** ERROR *** " fmt __RESET, __VA_ARGS__) 
#define error_puts(str)      fprintf(stderr, __BOLDRED "*** ERROR *** " str "\n" __RESET)

#define warning_printf(fmt, ...) fprintf(stdout, __YELLOW "*** WARNING *** " fmt __RESET, __VA_ARGS__)
#define warning_puts(str) puts(__YELLOW "*** WARNING *** " str __RESET)

#define ok_printf(fmt, ...) fprintf(stdout, __GREEN "--- OK --- " fmt __RESET, __VA_ARGS__)
#define ok_puts(str) puts(__GREEN "--- OK --- " str __RESET)

