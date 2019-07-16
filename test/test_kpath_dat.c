#include <kpath_dat.h>
#include <spg_wrap.h>
#include <poscar.h>
#include <error.h>

#define __DEBUG_ON
#include <debug.h>

#define PREC 1E-6

bool test_get_lattice_type(LatticeType latt_type);
void test_get_kpath();

int main(int argc, char* argv[]) {
  __UNUSED(argc);
  __UNUSED(argv);

  for (int i=1; i!=30; ++i) {
    test_get_lattice_type(i);
  }


  return 0;
}

bool test_get_lattice_type(LatticeType latt_type) {
  char filename[1 << 10] = {};
  char type_name[4] = {};

  kpt_print_lattice_type(latt_type, type_name);
  sprintf(filename, "../third_party/seekpath/seekpath/hpkot/band_path_data/%s/POSCAR_noinversion", type_name);

  POSCAR* poscar = psc_poscar_read(filename); 
  Cellp* cell = psc_poscar2spg_cell(poscar);

  const int hall_number = kpt_get_hall_number(cell, PREC);
  LatticeType analyzed_type = kpt_get_lattice_type(cell, PREC, hall_number);
  char analyzed_type_name[4];
  kpt_print_lattice_type(analyzed_type, analyzed_type_name);

  bool status = (latt_type == analyzed_type);

  if (status) {
    ok_printf("Testing get_lattice_type successed: expected %s, got %s\n", type_name, analyzed_type_name);
  } else {
    error_printf("Testing get_lattice_type failed: expected %s, but got %s\n", type_name, analyzed_type_name);
  }

  cellp_cell_free(cell);
  psc_poscar_free(poscar);

  return status;
}
