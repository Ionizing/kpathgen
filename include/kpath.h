#pragma once

#include <spglib/spglib.h>
#include "spg_wrap.h"

int kpt_get_lattice_type(const Cellp* cell,
                         const double prec,
                         const int    hall_number,
                               char   latt_type[4]);
