#define __DEBUG_ON
#include <debug.h>
#include <strops.h>

void test_str_n_tokens(const char* str) {
  debug_printf(" string = \"%s\"\n", str);
  debug_printf(" n_tokens = %d\n", str_n_tokens(str));
}

ERROR_CODE test_str_get_vecxi(char* str, int vec[], const int size) {
  ERROR_CODE status = str_get_vecxi(str, vec, size);
  debug_printf(" string = \"%s\"\n", str);
  debug_print_vecix(vec, size);
  return status;
}

ERROR_CODE test_str_get_vecxd(char* str, double vec[], const int size) {
  ERROR_CODE status = str_get_vecxd(str, vec, size); debug_printf(" string = \"%s\"\n", str);
  debug_print_vecdx(vec, size);
  return status;
}

bool test_str_start_with(const char* str, const char* pattern) {
  bool status = str_start_with(str, pattern);
  debug_printf(" string = \"%s\"\t pattern = \"%s\"\n\t\t\t judge result: %s\n", 
      str, pattern, true == status ? "true" : "false" );
  return status;
}

bool test_str_end_with(const char* str, const char* pattern) {
  bool status = str_end_with(str, pattern);
  debug_printf(" string = \"%s\"\t pattern = \"%s\"\n\t\t\t judge result: %s\n", 
      str, pattern, true == status ? "true" : "false" );
  return status;
}

void test_str_trim_left(const char* str, const char* seps) {
  char tmp[1 << 10];
  strcpy(tmp, str);
  debug_printf(" string = \"%s\"\t seps = \"%s\"\n", str, seps);
  char* result = str_trim_left(tmp, seps);
  debug_printf("\t result = \"%s\"\n", result);
  return;
}

void test_str_trim_right(const char* str, const char* seps) {
  char tmp[1 << 10];
  strcpy(tmp, str);
  debug_printf(" string = \"%s\"\t seps = \"%s\"\n", str, seps);
  char* result = str_trim_right(tmp, seps);
  debug_printf("\t result = \"%s\"\n", result);
  return;
}

void test_str_trim(const char* str, const char* seps) {
  char tmp[1 << 10];
  strcpy(tmp, str);
  debug_printf(" string = \"%s\"\t seps = \"%s\"\n", str, seps);
  char* result = str_trim(tmp, seps);
  debug_printf("\t result = \"%s\"\n", result);
  return;
}

void test_str_rev(const char* str) {
  char tmp[1 << 10];
  strcpy(tmp, str);
  str_rev(tmp);
  debug_printf(" string = \"%s\"\t reversed str = \"%s\"\n", str, tmp);
}


int main(int argc, char* argv[]) {
  __UNUSED(argc);
  __UNUSED(argv);


  warning_puts("test_str_n_tokens");
  test_str_n_tokens("abc def gl    sdfs w e ; sdfsljfsaoigoi3jr23 sdf; wr 23");
  test_str_n_tokens("asdf");
  test_str_n_tokens("");
  test_str_n_tokens("    )");
  test_str_n_tokens("\n13 \n 324 \n sdf");

  warning_puts("test_str_get_vecxi");
  char str[] = "1  2 3 4 5 6  7 8   9 0  ";
  int vec[10];
  test_str_get_vecxi(str, vec, 10);

  warning_puts("test_str_get_vecxd");
  double vec2[10];
  test_str_get_vecxd(str, vec2, 10);

  warning_puts("test_str_start_with");
  test_str_start_with("Foo Bar Baz", "Foo");
  test_str_start_with("Foo Bar Baz", " Foo");
  test_str_start_with("Foo Bar Baz", "");
  test_str_start_with("", "Foo");
  test_str_start_with("", "");

  warning_puts("test_str_end_with");
  test_str_end_with("Foo Bar Baz", "Baz");
  test_str_end_with("Foo Bar Baz", " Baz");
  test_str_end_with("Foo Bar Baz", "Foo");
  test_str_end_with("Foo", "Baz");
  test_str_end_with("", "");
  test_str_end_with("", "Baz");
  test_str_end_with("Baz", "Baz");
  test_str_end_with("Baz", "");

  warning_puts("test_str_trim_left");
  test_str_trim_left("   \tthis is a test string. .. \r\v", NULL);
  test_str_trim_right("   \tthis is a test string. .. \r\v", NULL);
  test_str_trim("   \tthis is a test string. .. \r\v", NULL);

  warning_puts("test_str_rev");
  test_str_rev("abc def ghi 123 456 789  ");
  test_str_rev("a");
  test_str_rev("");
  return 0;
}
