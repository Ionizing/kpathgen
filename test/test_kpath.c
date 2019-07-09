#include <kpath.h>
#include <poscar.h>
#include <error.h>
#include <spg_wrap.h>

#define __DEBUG_ON
#include <debug.h>

int main(int argc, char* argv[]) {
  __UNUSED(argc);
  __UNUSED(argv);

  POSCAR* poscar = psc_poscar_read("POSCAR");
  if (NULL == poscar) {
    error_puts("Parsing POSCAR failed.");
    return -1;
  } else {  }

  Cellp* cell = psc_poscar2spg_cell(poscar);
  if (NULL == cell) {
    error_puts("Converting from POSCAR to Cell failed.");
    return -1;
  } else {  }
  
  int res = kpt_standardize_cell(cell, true, false, 1E-5);
  if (0 == res) {
    error_puts("Cell standardization failed.");
    return -1;
  } else {  }

  psc_poscar_print(poscar);
  psc_poscar_free(poscar);
  poscar = spg_cell2psc_poscar(cell);
  psc_poscar_print(poscar);

  char symbol[7] = "";
  res = kpt_get_point_group(cell, 1E-5, symbol);
  if (0 == res) {
    error_puts("Get point group type failed.");
    return -1;
  } else {  } 

  ok_printf("The point group of current POSCAR is %s, point group number is %3d\n", symbol, res);


  psc_poscar_free(poscar);
  cellp_cell_free(cell);

  return 0;
}
