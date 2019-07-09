#pragma once

#include <ctype.h>

#include "base.h"
#include "error.h"
#include "debug.h"

/* count tokens separated by space, empty tokens are ignored */
int         str_n_tokens  (const char* str);

/* read double/int values from string and store them in a vector */
ERROR_CODE  str_get_vecxd (const char* str, double vec[], const int size);
ERROR_CODE  str_get_vecxi (const char* str, int    vec[], const int size);

/* check if string is start/end with given pattern */
bool        str_start_with(const char* str, const char* pattern);
bool        str_end_with  (const char* str, const char* patter);

/* trim string, return the trimmed string */
char*       str_trim_left (char* str, const char* seps);
char*       str_trim_right(char* str, const char* seps);
char*       str_trim      (char* str, const char* seps);

/* captalize/lower any alphabet characters in the string, then return it */
char*       str_toupper   (char* str);
char*       str_tolower   (char* str);

/* reverse a given string and return it */
char*       str_rev       (char* str);
