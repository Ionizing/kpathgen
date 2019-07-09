#pragma once

// standard libs
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// GNU Scientific Library libs
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

// spglib
#include "../third_party/spglib/src/spglib.h"

// Color definition
#define __RESET       "\033[0m"
#define __BLACK       "\033[30m"             /* Black */
#define __RED         "\033[31m"             /* Red */
#define __GREEN       "\033[32m"             /* Green */
#define __YELLOW      "\033[33m"             /* Yellow */
#define __BLUE        "\033[34m"             /* Blue */
#define __MAGENTA     "\033[35m"             /* Magenta */
#define __CYAN        "\033[36m"             /* Cyan */
#define __WHITE       "\033[37m"             /* White */
#define __BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define __BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define __BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define __BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define __BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define __BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define __BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define __BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


#define __STR_HELPER(x) #x
#define __STR(x) __STR_HELPER(x)

#define __VERSION_MAJOR 0
#define __VERSION_MINOR 0
#define __VERSION_PATCH 0
#define __VERSION "v" __STR(__VERSION_MAJOR) "." __STR(__VERSION_MINOR) "." __STR(__VERSION_PATCH) "\n"


#define __UNUSED(x) (void)x // disable unused parameter warning;
