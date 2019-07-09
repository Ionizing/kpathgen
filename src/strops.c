/* #define __DEBUG_ON */
#include "../include/debug.h"
#include "../include/strops.h"

int str_n_tokens(const char* str) {
  if (NULL == str || 0 == *str) {
    return 0;
  } else {  }

  debug_printf("Current string is: \"%s\"\n", str);

  int cnt = 0;
  char* to_be_dropped = (char*)calloc(1 << 15, sizeof(char));
  FILE* stream = fmemopen(str, strlen(str) + 1, "r");

  while (fscanf(stream, "%s", to_be_dropped) == 1) {
    ++cnt;
  }

  fclose(stream);
  free(to_be_dropped);
  return cnt;
}


ERROR_CODE str_get_vecxd(const char* str, double vec[], const int size) {
  if (size < 0) {
    error_puts("Invalid array size when reading from string.");
    return STRING_INVALID_VEC_SIZE; 
  } else {  }

  if (NULL == str || 0 == *str) {
    error_puts("Invalid string: NULL or EMPTY.");
    return STRING_NULL_OR_EMPTY;
  } else {  }

  debug_printf("Current string is: \"%s\"\n", str);
  /* linke string stream in cxx */
  FILE* stream = fmemopen(str, strlen(str) + 1, "r");

  for (int i=0; i!=size; ++i) {
    if (fscanf(stream, "%lf", &vec[i]) != 1) {
      error_puts("Invalid string to be parsed as float number.");
      /* debug_printf("current string is: \"%s\"", ) */
      return STRING_INVALID_STR_TO_FLOAT;
    } else {  }
  }

  fclose(stream);
  return SUCCESS;
}


ERROR_CODE str_get_vecxi(const char* str, int vec[], const int size) {
  if (size <= 0) {
    error_puts("Invalid array size when reading from string.");
    return STRING_INVALID_VEC_SIZE;
  } else if (!vec) {
    error_puts("Invalid vector pointer passed in: NULL pointer.");
    return STRING_NULL_OR_EMPTY;
  } else if (NULL == str || 0 == *str) {
    error_puts("Invalid string: NULL or EMPTY.");
    return STRING_NULL_OR_EMPTY;
  } else {  }
  
  debug_printf("Current string is: \"%s\"\n", str);

  FILE* stream = fmemopen(str, strlen(str) + 1, "r");

  for (int i=0; i!=size; ++i) {
    if (fscanf(stream, "%d", &vec[i]) != 1) {
      error_puts("Invalid string to be parsed as integer number.");
      return STRING_INVALID_STR_TO_INT;
    } else {  }
  }

  fclose(stream);
  return SUCCESS;
}


bool str_start_with(const char* str, const char* pattern) {
  if (!str || !pattern || !*str || !*pattern) {
    return false;
  }
  debug_printf("Current string is: \"%s\"\n", str);

  int lenA = strlen(str);
  int lenB = strlen(pattern);

  if (lenA < lenB) {
    return false;
  }

  const char* pA = str;
  const char* pB = pattern;
  while(*pA && *pB) {
    if (*pA != *pB) {
      return false;
    }
    ++pA;
    ++pB;
  }
  return true;
}

bool str_end_with(const char* str, const char* pattern) {
  if (!str || !pattern || !*str || !*pattern) {
    return false;
  } else {  }
  debug_printf("Current string is: \"%s\"\n", str);

  int lenA = strlen(str);
  int lenB = strlen(pattern);

  if (lenA < lenB) {
    return false;
  } else {  }

  const char* tailA = str + lenA - 1;
  const char* tailB = pattern + lenB - 1;

  while (tailA != str && tailB != pattern) {
    if (*tailA != *tailB) {
      return false;
    }
    --tailA;
    --tailB;
  }

  return true;
}

char* str_trim_left(char* str, const char* seps) {
  if (!str || !*str) {
    return NULL;
  } else {  }
  debug_printf("Current string is: \"%s\"\n", str);

  if (NULL == seps) {
    seps = "\t\n\v\f\r ";
  } else {  }

  size_t headlen = strspn(str, seps);
  if (headlen > 0) {
    size_t len = strlen(str);
    if (headlen == len) {
      str[0] = '\0';
    } else {
      memmove(str, str + headlen, len + 1 - headlen);
    }
  } else {  }

  return str;
}


char* str_trim_right(char* str, const char* seps) {
  if (!str || !*str) {
    return str;
  } else {  }
  debug_printf("Current string is: \"%s\"\n", str);

  if (NULL == seps) {
    seps = "\t\n\v\f\r ";
  } else { }

  char* ptail = str + strlen(str) - 1;
  while (str != ptail && NULL != strchr(seps, *ptail)) {
    *ptail = '\0';
    --ptail;
  }
  return str;
}


inline char* str_trim(char* str, const char* seps) {
  return str_trim_left(str_trim_right(str, seps), seps);
}


char* str_toupper(char* str) {
  if (!str || !*str) {
    return str;
  } else {  }
  debug_printf("Current string is: \"%s\"\n", str);

  char* ptr = str;   // backup head pointer
  while ('\0' != *ptr) {
    *ptr = toupper(*ptr);
    ++ptr;
  }

  return str;
}


char* str_tolower(char* str) {
  if (!str || !*str) {
    return str;
  } else {  }
  debug_printf("Current string is: \"%s\"\n", str);

  char* ptr = str;   // backup head pointer
  while ('\0' != *ptr) {
    *ptr = tolower(*ptr);
    --ptr;
  }

  return str;
}


char* str_rev(char* str) {
  if (NULL == str || '\0' == *str) {
    return str;
  }
  debug_printf("Current string is: \"%s\"\n", str);

  char* tail = str;
  char* head = str;
  while (*tail) { ++tail; }  // now tail is pointing the '\0'
  --tail;
  char tmpch;
  while (head < tail) {
    tmpch = *head;
    *head = *tail;
    *tail = tmpch;
    ++head;
    --tail;
  }

  return str;
}
