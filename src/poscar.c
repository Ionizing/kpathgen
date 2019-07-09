/* #define __DEBUG_ON */
#include "../include/poscar.h"
#include "../include/debug.h"

POSCAR* psc_poscar_alloc(const int n_atoms, const int n_types) {
  
  if (n_atoms < 1 || n_types < 1) {
    error_puts("Bad poscar alloc request: n_atoms or n_types < 1.");
    return NULL;
  } else {  }

  POSCAR* poscar = NULL;
  if ((poscar = (POSCAR*)calloc(1, sizeof(POSCAR))) == NULL) {
    error_puts("Memory could not be allocated for poscar.");
    poscar = NULL;
    return NULL;
  }

  poscar->n_atoms = n_atoms;
  poscar->n_types = n_types;

  if ((poscar->types = (int*)calloc(n_types, sizeof(int))) == NULL) {
    error_puts("Memory could not be allocated for poscar.");
    poscar->types = NULL;
    free(poscar);
    poscar = NULL;
    return NULL;
  } else {  }

  if ((poscar->type_names = (char (*)[10])calloc(n_types, sizeof(char[10]))) == NULL) {
    error_puts("Memory could not be allocated for poscar.");
    poscar->type_names = NULL;
    free(poscar->types);
    poscar->types = NULL;
    free(poscar);
    poscar = NULL;
    return NULL;
  } else {  }


  if ((poscar->positions = 
        (double (*)[3])calloc(n_atoms, sizeof(double[3]))) == NULL) {
    error_puts("Memory could not be allocated for poscar.");
    free(poscar->types);
    poscar->types = NULL;
    free(poscar->type_names);
    poscar->type_names = NULL;
    free(poscar);
    poscar = NULL;
    return NULL;
  } else {  }

  *(poscar->comment) = '\0';

  return poscar;
} // end of psc_poscar_alloc


ERROR_CODE psc_poscar_free(POSCAR* poscar) {
  if (NULL != poscar) {
    if (NULL != poscar->positions) {
      free(poscar->positions);
      poscar->positions = NULL;
    }
    if (NULL != poscar->types) {
      free(poscar->types);
      poscar->types = NULL;
    }
    if (NULL != poscar->type_names) {
      free(poscar->type_names);
      poscar->type_names = NULL;
    }
    free(poscar);
  }
  return SUCCESS;
} // end of psc_poscar_free


ERROR_CODE psc_poscar_set(POSCAR* poscar,
                          const char* comment,
                          const double lattice[3][3],
                          const int types[],
                          const char type_names[][10],
                          const int n_types,
                          const double positions[][3],
                          const int n_atoms) {
  if (NULL == poscar) {
    error_puts("NULL poscar object pointer passed in when setting it.");
    return POSCAR_NULL_OBJECT;
  } else {  }

  if (n_types != poscar->n_types || n_atoms != poscar->n_atoms) {
    error_puts("Inconsistent n_types or n_atoms with given poscar object.");
    return POSCAR_INVALID_ATOM_NUMBER; 
  } else {  }

  memcpy(poscar->lattice, lattice, 9 * sizeof(double));
  debug_print_matd33(lattice);
  debug_print_matd33(poscar->lattice);

  strncpy(poscar->comment, comment, sizeof(poscar->comment));

  debug_print_vecstr(type_names, n_types);
  memcpy(poscar->type_names, type_names, sizeof(char[10]) * n_types);
  debug_print_vecstr(poscar->type_names, n_types);

  debug_print_vecix(types, n_types);
  memcpy(poscar->types, types, n_types * sizeof(int));
  debug_print_vecix(poscar->types, n_types);

  memcpy(poscar->positions, positions, 3 * n_atoms * sizeof(double));
  debug_print_matdx3(positions, n_atoms);
  debug_print_matdx3(poscar->positions, n_atoms);

  return SUCCESS;
} // end of psc_poscar_set


POSCAR* psc_poscar_read(const char* filename) {
  /* prepare the buffer  */
  static const int bufsize = 5 << 20;

  /* prepare the temporary poscar data */ 
  char comment[1 << 10] = {0};
  double lattice[3][3];
  int n_types = 0;
  int* types;
  char (*type_names)[10];
  int n_atoms = 0;
  double (*positions)[3];



  char* buf = (char*)calloc(bufsize, sizeof(char)); // allocate 1 MB buffer
  if (NULL == buf) {
    error_puts("Memory cannot be allocated for poscar.");
    return NULL;
  }

  /* manually specify the buffer */
  setvbuf(stdin, buf, _IOFBF, bufsize);

  /* open file */
  FILE* fp = fopen(filename, "r");
  if (NULL == fp) {
    error_printf("Open POSCAR file \"%s\" failed!\n", filename);
    return NULL;
  }

  ok_printf("Opened POSCAR file \"%s\" successfully.\n", filename);

  char* line     = NULL;  // getline will reallocate memory for `line`
  size_t linelen = 0;
  ssize_t nread;        // number of read bytes

  /* start parsing */
  if((nread = getline(&line, &linelen, fp)) != -1) {
    str_trim_right(line, NULL);
    strncpy(comment, line, (1 << 10) - 1); // comment line
  } else {
    /* error_puts("Empty file.") */
  }
  
  /* scale constant */
  double scale = 1.0;
  if (1 != fscanf(fp, "%lf", &scale)) {
    debug_puts("here");
    error_puts("Cannot read a float variable from POSCAR file.");
    return NULL;
  } else if (scale <= 1E-5) {
    debug_puts("here");
    error_puts("Invalid lattice scale parameter in POSCAR.");
    return NULL;
  } else {  }
  
  getline(&line, &linelen, fp); // skip the newline character
  
  /* basic vectors */
  for (int i=0; i!=3; ++i) {
    if ((nread = getline(&line, &linelen, fp)) != -1) {
      str_get_vecxd(line, lattice[i], 3);
    } else {
      /* error_puts  */
    }
  }

  mat_mat33_multiply(lattice, scale);

  /* read element tag, this part is indispensable and cannot be missing */
  if ((nread = getline(&line, &linelen, fp)) != -1) {
    str_trim(line, NULL);
  } else { 
    /* error_puts("Incomplete content in POSCAR.") */
  }

  n_types = str_n_tokens(line); // atom types.
  type_names = (char (*)[10]) calloc(n_types, sizeof(char[10]));

  { // Tricky method here, surround this block with `{}` to release variables automaticaly
    FILE* type_str_fp = fmemopen(line, 1 << 10, "r");
    for (int i=0; i!=n_types; ++i) {
      fscanf(type_str_fp, "%9s", type_names[i]);
      debug_printf("type names = \"%s\"\n", type_names[i]);
      if (!isalpha(type_names[i][0])) {
        error_printf("Invalid element tag: \"%s\"\n", type_names[i]);
        return NULL;
      } else if (strlen(type_names[i]) >= 9) {
        error_puts("Element tag length too long, less than 9 is allowed.");
        return NULL;
      } 
    }
    fclose(type_str_fp);
  }

  /* start parsing the number of atoms correspoing to type_names */
  if ((nread = getline(&line, &linelen, fp)) != -1) {
    debug_printf("number of atoms is: \"%s\"\n", line);
    str_trim(line, NULL);
  } else {
    /* error_puts("Incomplete content in POSCAR.")  */
  }

  if (!line || !isnumber(line[0])) { // Line contains no numbers
    error_printf("No number of atoms corresponding to atom type tags. \n\tCurrent line: \"%s\"\n", line);
    return NULL;
  } else if (n_types != str_n_tokens(line)) { // atom quantities
    error_printf("Atom types' number inconsistent with atom type quantities number. \n\tCurrent line: \"%s\"\n", line);
  } else {  }
  
  types = (int*) calloc(n_types, sizeof(int));
  str_get_vecxi(line, types, n_types);
  n_atoms = 0;
  for (int i=0; i!=n_types; ++i) {
    n_atoms += types[i];      // get number of atoms
  }
  
  if ((nread = getline(&line, &linelen, fp)) != -1) {
    str_trim(line, NULL);
  } else {
    error_puts("Incomplete content in POSCAR.");
    debug_puts("here");
    return NULL;
  }

  if (!line || !isalpha(line[0])) {
    error_printf("Invalid content: cartesian or direct tag needed.\n\tCurrent line: \"%s\"\n", line);
    return NULL;
  } else if ('S' == toupper(line[0])) {
    if ((nread = getline(&line, &linelen, fp)) != -1) {
      str_trim(line, NULL);
    } else {
      error_puts("Incomplete content in POSCAR.");
      return NULL;
    }
  }

  bool is_frac_coord = false;

  switch (line[0]) {
    case 'c': ;
    case 'C': is_frac_coord = false;
              break;

    case 'd': ;
    case 'D': is_frac_coord = true;
              break;

    default: error_printf("Incomplete content in POSCAR: No coordinate type tag.\n\tCurrent line: \"%s\"", line);
             return NULL;
  }

  /* start reading atom positions, waited too long >_< */
  positions = (double (*)[3]) calloc(n_atoms, sizeof(double [3]));
  double bcell[3][3];
  mat_mat33_inverse(lattice, bcell);

  for (int i=0; i!=n_atoms; ++i) {
    getline(&line, &linelen, fp);
    double tmppos[3];
    str_get_vecxd(line, tmppos, 3);
    if (false == is_frac_coord) {
      mat_bra3_dot_mat33(tmppos, bcell, positions[i]);
    } else {
      memcpy(positions[i], tmppos, sizeof(tmppos));
    }
  }

  /*  convert row major to col major  */
  /* No need. */
  /* mat_mat33_trans(lattice); */
  
  /* Now move the parsed data to poscar object */
  POSCAR* poscar = psc_poscar_alloc(n_atoms, n_types);

  psc_poscar_set(poscar, comment, lattice, types, type_names, n_types, positions, n_atoms);

  free(types);
  free(positions);

  /*
   * when repeatly invoking `getline`, it calls `realloc`, which frees the old block of memory
   * that `ptr` pointed to, if needed. However, user should manually free `ptr` if no more
   * `getline` is invoked.
   */

  free(line);
  free(buf);
  fclose(fp);
  return poscar;
} // end of psc_poscar_read


/* Auxiliary function: write poscar content to a file pointer */
void _psc_poscar_write_to_fp(const POSCAR* poscar, FILE* fp) {
  /* I will not check if `fo` is open in "w" mode, because this function will not
   * be in public */
  fputs(poscar->comment, fp); fputs("\n", fp);
  fputs("\t1.00000\n", fp);
  fprintf(fp, " %10.5f %10.5f %10.5f\n", 
      poscar->lattice[0][0], poscar->lattice[0][1], poscar->lattice[0][2]);
  fprintf(fp, " %10.5f %10.5f %10.5f\n", 
      poscar->lattice[1][0], poscar->lattice[1][1], poscar->lattice[1][2]);
  fprintf(fp, " %10.5f %10.5f %10.5f\n", 
      poscar->lattice[2][0], poscar->lattice[2][1], poscar->lattice[2][2]);

  for (int i=0; i!=poscar->n_types; ++i) {
    fprintf(fp, "%5s", poscar->type_names[i]);
  } fputs("\n", fp);

  for (int i=0; i!=poscar->n_types; ++i) {
    fprintf(fp, "%5d", poscar->types[i]);
  } fputs("\n", fp);

  /* Ignore the selective dynamics tag */
  fputs("Direct\n", fp);

  for (int i=0; i!=poscar->n_atoms; ++i) {
    fprintf(fp, "  %10.6f %10.6f %10.6f\n", 
        poscar->positions[i][0], poscar->positions[i][1], poscar->positions[i][2]);
  }

}
/* Auxiliary function ends here */


ERROR_CODE psc_poscar_write(const POSCAR* poscar,
                            const char* filename) {
  if (NULL == poscar) {
    error_puts("Invalid poscar object passed in: NULL object.");
    return POSCAR_NULL_OBJECT;
  } else {  }

  FILE* fp = fopen(filename, "w");
  if (NULL == fp) {
    error_printf("Cannot write to \"%s\": open failed,\n", filename);
    return FILE_OPEN_FAILED;
  }

  _psc_poscar_write_to_fp(poscar, fp);

  fclose(fp);
  return SUCCESS;
} // end of psc_poscar_write


ERROR_CODE psc_poscar_print(const POSCAR* poscar) {
  if (NULL == poscar) {
    error_puts("Invalid poscar object passed in: NULL object.");
    return POSCAR_NULL_OBJECT;
  } else {  }

  fprintf(stdout, __CYAN "This is the content of current poscar object.\n" __RESET);
  fprintf(stdout, __GREEN);
  _psc_poscar_write_to_fp(poscar, stdout);
  fprintf(stdout, __RESET);

  return SUCCESS;
} // end of psc_poscar_print
