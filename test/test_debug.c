#define __DEBUG_ON
#include <debug.h>
#include <error.h>

int main(int argc, char* argv[]) {
  __UNUSED(argc);
  __UNUSED(argv);

  error_puts("test error puts.");
  error_printf("test error %d\n", 1);

  warning_puts("test warning puts.");
  warning_printf("test warning %d\n", 2);

  double mat_d3[3][3] = {
    {1, 2, 3},
    {4, 5, 6},
    {7, 8, 8} };

  debug_print_matd33(mat_d3);

  int mat_i3[3][3] = {
    {1, 2, 3},
    {4, 5, 6},
    {7, 8, 9} };

  debug_print_mati33(mat_i3);

  double vec_d3[3] = {0, 1, 2};

  debug_print_vecd3(vec_d3);

  int vec_i3[3] = {0, 1, 2};

  debug_print_veci3(vec_i3);

  double positions[9][3];
  for (int i=0; i!=9; ++i) {
    for (int j=0; j!=3; ++j) {
      positions[i][j] = i*10 + j;
    }
  }
  debug_print_matdx3(positions, 9);
  debug_print_vecdx(vec_d3, 3);
  debug_print_vecix(vec_i3, 3);

  const char vecstr[][10] = {"a", "abc", "str", "114514", "810", "8*7", "8*9"};
  warning_puts("test_debug_vecstr");
  debug_print_vecstr(vecstr, 7);
  return 0;
}
