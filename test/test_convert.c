#define __DEBUG_ON

#include <spg_wrap.h>
#include <debug.h>

int main(int argc, char* argv[]) {
  __UNUSED(argc);
  __UNUSED(argv);

  POSCAR* poscar = psc_poscar_read("POSCAR");
  if (NULL == poscar) {
    error_puts("Parsing POSCAR failed.");
    return -1;
  } else {  }

  psc_poscar_print(poscar);

  Cellp* cell = psc_poscar2spg_cell(poscar);
  if (NULL == cell) {
    error_puts("Converting POSCAR to Cell failed");
    return -1;
  } else {  }

  debug_print_matd33(cell->lattice); 
  debug_print_vecix(cell->types, cell->size);
  debug_print_matdx3(cell->position, cell->size);

  psc_poscar_free(poscar);
  poscar = spg_cell2psc_poscar(cell);
  if (NULL == poscar) {
    error_puts("Converting Cell to POSCAR failed.");
    return -1;
  } else { }

  psc_poscar_print(poscar);
  psc_poscar_free(poscar);
  cellp_cell_free(cell);

  return 0;
}
