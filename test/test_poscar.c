#define __DEBUG_ON

#include <poscar.h>
#include <debug.h>

int main(int argc, char* argv[]) {
  __UNUSED(argc);
  __UNUSED(argv);

  const int n_atoms = 10;
  const int n_types = 2;
  POSCAR* poscar = psc_poscar_alloc(n_atoms, n_types);
  /* debug_printf("%s:%d:%s %s = %p\n", __FILE__, __LINE__, __func__, __STR(poscar), poscar); */
  const char* comment = "This is a comment";
  const double lattice[3][3] = {
    {1, 2, 3},
    {0, 1, 0},
    {0, 0, 1} };
  const int types[] = {4, 6};
  const char type_names[][10] = {"C", "H"};
  /* debug_print_vecix(types, 10); */

  const double positions[][3] = {
    {0, 0, 0},
    {0, 0, 1},
    {0, 1, 0},
    {0, 1, 1},
    {1, 0, 0},
    {1, 0, 1},
    {1, 1, 0},
    {1, 1, 1},
    {1, 1, 1},
    {0, 0, 0}
  };
  /* debug_print_matdx3(positions, 10); */
  
  psc_poscar_set(poscar, comment, lattice, types, type_names, n_types, positions, n_atoms);
  psc_poscar_print(poscar);
  psc_poscar_free(poscar);
  
  if (NULL == (poscar = psc_poscar_read("./POSCAR"))) {
    error_puts("Open POSCAR failed.");
  } else {
    psc_poscar_print(poscar);
    psc_poscar_write(poscar, "POSCAR_parsed.vasp");
    psc_poscar_free(poscar);
  }

  return 0;
}
