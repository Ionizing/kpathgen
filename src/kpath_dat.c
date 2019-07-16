#include "../include/kpath_dat.h"

/* Auxiliary Functions */
typedef struct {
  double  lattice[3][3];  // ROW MAJOR
  double rlattice[3][3];  // basis of lattice in reciprocal space
  double  a,  b,  c;
  double ra, rb, rc;      // lengths of reciprocal basis
  double  alpha,  beta,  gamma; // angels of real basis
  double ralpha, rbeta, rgamma; // angels of reciprocal basis
} _lattice_detail;

typedef SpacegroupType spg_SpacegroupType;


void _calc_lattice_detail(_lattice_detail* _l) {
  if (NULL == _l) {
    return ;
  }

  _l->a = sqrt(
       (_l->lattice[0][0]) * (_l->lattice[0][0]) +
       (_l->lattice[0][1]) * (_l->lattice[0][1]) +
       (_l->lattice[0][2]) * (_l->lattice[0][2]) );

  _l->b = sqrt(
       (_l->lattice[1][0]) * (_l->lattice[1][0]) +
       (_l->lattice[1][1]) * (_l->lattice[1][1]) +
       (_l->lattice[1][2]) * (_l->lattice[1][2]) );

  _l->c = sqrt(
       (_l->lattice[2][0]) * (_l->lattice[2][0]) +
       (_l->lattice[2][1]) * (_l->lattice[2][1]) +
       (_l->lattice[2][2]) * (_l->lattice[2][2]) );

  mat_mat33_inverse(_l->lattice, _l->rlattice);
  _l->ra = sqrt(
       (_l->rlattice[0][0]) * (_l->rlattice[0][0]) +
       (_l->rlattice[0][1]) * (_l->rlattice[0][1]) +
       (_l->rlattice[0][2]) * (_l->rlattice[0][2]) );
  
  _l->rb = sqrt(
       (_l->rlattice[1][0]) * (_l->rlattice[1][0]) +
       (_l->rlattice[1][1]) * (_l->rlattice[1][1]) +
       (_l->rlattice[1][2]) * (_l->rlattice[1][2]) );
  
  _l->rc = sqrt(
       (_l->rlattice[2][0]) * (_l->rlattice[2][0]) +
       (_l->rlattice[2][1]) * (_l->rlattice[2][1]) +
       (_l->rlattice[2][2]) * (_l->rlattice[2][2]) );
  
  _l->alpha = acos((_l->lattice[1][0] * _l->lattice[2][0] + 
                    _l->lattice[1][1] * _l->lattice[2][1] +
                    _l->lattice[1][2] * _l->lattice[2][2]) ) / (_l->b * _l->c);

  _l->beta  = acos((_l->lattice[1][0] * _l->lattice[2][0] + 
                    _l->lattice[1][1] * _l->lattice[2][1] +
                    _l->lattice[1][2] * _l->lattice[2][2]) ) / (_l->a * _l->c);

  _l->gamma = acos((_l->lattice[1][0] * _l->lattice[2][0] + 
                    _l->lattice[1][1] * _l->lattice[2][1] +
                    _l->lattice[1][2] * _l->lattice[2][2]) ) / (_l->a * _l->b);

  _l->ralpha = acos((_l->rlattice[1][0] * _l->rlattice[2][0] + 
                     _l->rlattice[1][1] * _l->rlattice[2][1] +
                     _l->rlattice[1][2] * _l->rlattice[2][2]) ) / (_l->rb * _l->rc);

  _l->rbeta  = acos((_l->rlattice[1][0] * _l->rlattice[2][0] + 
                     _l->rlattice[1][1] * _l->rlattice[2][1] +
                     _l->rlattice[1][2] * _l->rlattice[2][2]) ) / (_l->ra * _l->rc);

  _l->rgamma = acos((_l->rlattice[1][0] * _l->rlattice[2][0] + 
                     _l->rlattice[1][1] * _l->rlattice[2][1] +
                     _l->rlattice[1][2] * _l->rlattice[2][2]) ) / (_l->ra * _l->rb);

  

}


/* Auxiliary Functions ends here */


LatticeType kpt_get_lattice_type(const Cellp* cell,
                                 const double prec,
                                 const int    hall_number) {

  typedef enum {
    _placeholder = 0,
    _cP, _cF, _cI, /**/ /**/ /**/ /**/
    _tP, /**/ _tI, /**/ /**/ /**/ /**/
    _oP, _oF, _oI, /**/ _oC, _oA, /**/
    _hP, /**/ /**/ /**/ /**/ /**/ _hR,
    _mP, /**/ /**/ /**/ _mC, /**/ /**/
    _aP, /**/ /**/ /**/ /**/ /**/ /**/ 
  } _LatticeType_Coarse;

  static const _LatticeType_Coarse lattice_type_table[7][9] = {
    {   0,   0,   0,   0,   0,   0,   0,   0},
    {   0, _cP, _cF, _cI,   0,   0,   0,   0},
    {   0, _tP,   0, _tI,   0,   0,   0,   0},
    {   0, _oP, _oF, _oI,   0, _oC, _oA,   0},
    {   0, _hP,   0,   0,   0,   0,   0, _hR},
    {   0, _mP,   0,   0,   0, _mC,   0,   0},
    {   0, _aP,   0,   0,   0,   0,   0,   0},
  };


  spg_SpacegroupType spacegroup = spgdb_get_spacegroup_type(hall_number);
  const int spacegroup_number = spacegroup.number;
  const int centering = (int)spacegroup.centering;  // cast enum Centering into int

  int lattice_system = 0;
  switch (spacegroup_number) {
    case    1 ...   2 : lattice_system = 6; break;
    case    3 ...  15 : lattice_system = 5; break;
    case   16 ...  74 : lattice_system = 3; break;
    case   75 ... 142 : lattice_system = 2; break;
    case  143 ... 194 : lattice_system = 4; break;
    case  195 ... 230 : lattice_system = 1; break;
  }

  const _LatticeType_Coarse lattice_type_coarse = lattice_type_table[lattice_system][centering];

  switch (lattice_type_coarse) {
    case _cI: return _cI1;
    case _tP: return _tP1;
    case _oP: return _oP1;
    case _mP: return _mP1;
    case _cP: return (spacegroup_number <= 206) ? _cP1 : _cP2;
    case _cF: return (spacegroup_number <= 206) ? _cF1 : _cF2;

    default: break;
  } 

  _lattice_detail _l;
  memcpy(_l.lattice, cell->lattice, sizeof(double[3][3]));
  mat_mat33_trans(_l.lattice);
  _calc_lattice_detail(&_l);

  switch (lattice_type_coarse) {
    case _tI: return (_l.c < _l.a) ? _tI1 : _tI2;
    case _oF: {
      const double ad2 = 1.0 / (_l.a * _l.a);
      const double bd2 = 1.0 / (_l.b * _l.b);
      const double cd2 = 1.0 / (_l.c * _l.c);
      return (ad2 > bd2 + cd2) ? _oF1 :               // 1/a2 > 1/b2 + 1/c2
            ((cd2 > ad2 + bd2) ? _oF2 :               // 1/c2 > 1/a2 + 1/b2
                                 _oF3);               // 1/b2 > 1/a2 + 1/c2
    }
    case _oI: {
      return (_l.c > _l.a && _l.c > _l.b) ? _oI1 :    // c largest
            ((_l.b > _l.a && _l.b > _l.c) ? _oI2 :    // b largest
                                            _oI3);    // a largest
    }
    case _oC: return (_l.a < _l.b) ? _oC1 : _oC2;
    case _oA: return (_l.b < _l.c) ? _oA1 : _oA2;
    case _hP: {
      switch (spacegroup_number) {
        case 143 ... 149:
        case 151:
        case 153:
        case 157:
        case 159 ... 163:  return _hP1;
        default:           return _hP2;
      }         
    }
    case _hR: return (_l.a < 0.8164965809 * _l.c) ? _hR1 : _hR2;  // sqrt(3) * a < sqrt(2) * c ?
    case _mC: return (_l.b < _l.a * sin(_l.beta)) ? _mC1 :
              ( ( (_l.a*_l.a * sin(_l.beta)*sin(_l.beta) / (_l.b*_l.b) - 
                   _l.a*cos(_l.beta)/_l.c) < 1 )? _mC2 : _mC3);
    case _aP: return (_l.ralpha < 0 && _l.rbeta < 0 && _l.rgamma < 0) ? _aP2 : _aP3;
              // interaxial angels are acute -> _aP2, else _aP3
    default: return _null;
  }

  return _null;
} // end of kpt_get_lattice_type


ERROR_CODE kpt_print_lattice_type(const LatticeType latt_type,
                                        char        res_str[4]) {
  static const char type_strs[][4] = {
    "NUL", "cP1", "cP2", "cF1",
    "cF2", "cI1", "tP1", "tI1",
    "tI2", "oP1", "oF1", "oF2",
    "oF3", "oI1", "oI2", "oI3",
    "oC1", "oC2", "oA1", "oA2",
    "hP1", "hP2", "hR1", "hR2",
    "mP1", "mC1", "mC2", "mC3",
    "aP2", "aP3",
  };

  if (NULL == res_str) {
    error_puts("NULL string pointer passed in when printing lattice type.");
    return STRING_NULL_OR_EMPTY;
  } else {  }

  strcpy(res_str, type_strs[latt_type]);
  return SUCCESS;
}


ERROR_CODE kpt_print_high_symmetry_point(const HighSymmetryPoint hsp,
                                               char              res_str[9]) {
  static const char hsp_strs[][9] = {
         "A",      "A_0",
         "B",      "B_0",      "B_2",
         "C",      "C_0",      "C_2",      "C_4",
         "D",      "D_0",      "D_2",      
               "DELTA_0",                          
         "E",      "E_0",      "E_2",      "E_4",
         "F",      "F_0",      "F_2",      "F_4",
         "G",      "G_0",      "G_2",      "G_4",      "G_6",
     "GAMMA",                                    
         "H",      "H_0",      "H_2",      "H_4",      "H_6",
         "I",                  "I_2",                          
                   "J_0",                                    
         "K",                  "K_2",      "K_4",                
         "L",      "L_0",      "L_2",      "L_4",      
  "LAMBDA_0",                                    
         "M",      "M_0",      "M_2",      "M_4",      "M_6",     "M_8",
         "N",      "N_2",      "N_4",      "N_6",      
         "P",      "P_0",      "P_2",                
                   "Q_0",                                    
         "R",      "R_0",      "R_2",                
         "S",      "S_0",      "S_2",      "S_4",      "S_6",
               "SIGMA_0",                          
         "T",                  "T_2",                
         "U",      "U_0",      "U_2",      
         "V",      "V_0",      "V_2",      
         "W",                  "W_2",                
         "X",                  "X_1",                
         "Y",      "Y_0",      "Y_2",      "Y_4",
         "Z",      "Z_0",      "Z_2",
  };

  if (NULL == res_str) {
    error_puts("NULL string pointer passed in when print high symmetry point.");
    return STRING_NULL_OR_EMPTY;
  }

  strcpy(res_str, hsp_strs[hsp]);
  return SUCCESS;
}


/*
 * KPath data
 */
static const KPath kpath_dat[30] = { // Only takes up 65 kB
  [_null] = { // null type
    .latt_type        = _null,
    .is_need_calculated = false,
    .n_highsym_points = 0,
    .n_paths          = 0,
    .highsym_points   = {},
    .path             = {},
    .hsp_coordinates  = {},
  },
  [_cP1] = { // cP1
    .latt_type        = _cP1,
    .is_need_calculated = false,
    .n_highsym_points = 5,
    .n_paths          = 7,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_R,
      _hsp_M,
      _hsp_X,
      _hsp_X_1,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X    },
      {_hsp_X,        _hsp_M    },
      {_hsp_M,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_R    },
      {_hsp_R,        _hsp_X    },
      {_hsp_R,        _hsp_M    },
      {_hsp_M,        _hsp_X_1  },
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_R]      = {   0.5,   0.5,   0.5},
      [_hsp_M]      = {   0.5,   0.5,     0},
      [_hsp_X]      = {     0,   0.5,     0},
      [_hsp_X_1]    = {   0.5,     0,     0},
    },
  },
  [_cP2] = { // cP2
    .latt_type        = _cP2,
    .is_need_calculated = false,
    .n_highsym_points = 5,
    .n_paths          = 6,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_R,
      _hsp_M,
      _hsp_X,
      _hsp_X_1,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X    },
      {_hsp_X,        _hsp_M    },
      {_hsp_M,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_R    },
      {_hsp_R,        _hsp_X    },
      {_hsp_R,        _hsp_M    },
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_R]      = {   0.5,   0.5,   0.5},
      [_hsp_M]      = {   0.5,   0.5,     0},
      [_hsp_X]      = {     0,   0.5,     0},
      [_hsp_X_1]    = {   0.5,     0,     0},
    },
  },
  [_cF1] = { // cF1
    .latt_type        = _cF1,
    .is_need_calculated = false,
    .n_highsym_points = 7,
    .n_paths          = 7,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_X,
      _hsp_L,
      _hsp_W,
      _hsp_W_2,
      _hsp_K,
      _hsp_U,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_U},
      {_hsp_K,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_L},
      {_hsp_L,        _hsp_W},
      {_hsp_W,        _hsp_X},
      {_hsp_X,        _hsp_W_2}, 
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_X]      = {   0.5,     0,   0.5},
      [_hsp_L]      = {   0.5,   0.5,   0.5},
      [_hsp_W]      = {   0.5,  0.25,  0.75},
      [_hsp_W_2]    = {  0.75,  0.25,   0.5},
      [_hsp_K]      = { 0.375, 0.375,  0.75},
      [_hsp_U]      = { 0.635,  0.25, 0.625},
    },
  },
  [_cF2] = { // cF2
    .latt_type        = _cF2,
    .is_need_calculated = false,
    .n_highsym_points = 7,
    .n_paths          = 6,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_X,
      _hsp_L,
      _hsp_W,
      _hsp_W_2,
      _hsp_K,
      _hsp_U,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_U},
      {_hsp_K,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_L},
      {_hsp_L,        _hsp_W},
      {_hsp_W,        _hsp_X},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_X]      = {   0.5,     0,   0.5},
      [_hsp_L]      = {   0.5,   0.5,   0.5},
      [_hsp_W]      = {   0.5,  0.25,  0.75},
      [_hsp_W_2]    = {  0.75,  0.25,   0.5},
      [_hsp_K]      = { 0.375, 0.375,  0.75},
      [_hsp_U]      = { 0.625,  0.25, 0.625},
    },
  },
  [_cI1] = { // cI1
    .latt_type        = _cI1,
    .is_need_calculated = false,
    .n_highsym_points = 4,
    .n_paths          = 6,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_H,
      _hsp_P,
      _hsp_N,
    },
    .path = {
      {_hsp_GAMMA,    _hsp_H},
      {_hsp_H,        _hsp_N},
      {_hsp_N,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_P},
      {_hsp_P,        _hsp_H},
      {_hsp_P,        _hsp_N},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_H]      = {   0.5,   -.5,   0.5},
      [_hsp_P]      = {  0.25,  0.25,  0.25},
      [_hsp_N]      = {     0,     0,   0.5},
    },
  },
  [_tP1] = { // tP1
    .latt_type        = _tP1,
    .is_need_calculated = false,
    .n_highsym_points = 6,
    .n_paths          = 9,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_Z,
      _hsp_M,
      _hsp_A,
      _hsp_R,
      _hsp_X,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_M},
      {_hsp_M,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_Z,        _hsp_R},
      {_hsp_R,        _hsp_A},
      {_hsp_A,        _hsp_Z},
      {_hsp_X,        _hsp_R},
      {_hsp_M,        _hsp_A},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_Z]      = {     0,     0,   0.5},
      [_hsp_M]      = {   0.5,   0.5,     0},
      [_hsp_A]      = {   0.5,   0.5,   0.5},
      [_hsp_R]      = {     0,   0.5,   0.5},
      [_hsp_X]      = {     0,   0.5,     0},
    },
  },
  [_tI1] = { // tI1  --> Need to be calculated
    .latt_type        = _tI1,
    .is_need_calculated = true,
    .n_highsym_points = 7,
    .n_paths          = 8,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_M,
      _hsp_X,
      _hsp_P,
      _hsp_Z,
      _hsp_Z_0,
      _hsp_N,
    },
    .path = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_M},
      {_hsp_M,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_Z_0,      _hsp_M},
      {_hsp_X,        _hsp_P},
      {_hsp_P,        _hsp_N},
      {_hsp_N,        _hsp_GAMMA},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_M]      = {   -.5,   0.5,   0.5},
      [_hsp_X]      = {     0,     0,   0.5},
      [_hsp_P]      = {  0.25,  0.25,  0.25},
      [_hsp_Z]      = {    -1,    -1,    -1},  // Need to be calculated.
      [_hsp_Z_0]    = {    -1,    -1,    -1},  // Need to be calculated.
      [_hsp_N]      = {     0,   0.5,     0},
    },
  },
  [_tI2] = { // tI2 --> Need to be calculated
    .latt_type        = _tI2,
    .is_need_calculated = true,
    .n_highsym_points = 9,
    .n_paths          = 9,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_M,
      _hsp_X,
      _hsp_P,
      _hsp_N,
      _hsp_S_0,
      _hsp_S,
      _hsp_R,
      _hsp_G,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_P},
      {_hsp_P,        _hsp_N},
      {_hsp_N,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_M},
      {_hsp_M,        _hsp_S},
      {_hsp_S_0,      _hsp_GAMMA},
      {_hsp_X,        _hsp_R},
      {_hsp_G,        _hsp_M},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_M]      = {   0.5,   0.5,   0.5},
      [_hsp_X]      = {     0,     0,   0.5},
      [_hsp_P]      = {  0.25,  0.25,  0.25},
      [_hsp_N]      = {     0,   0.5,     0},
      [_hsp_S_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_S]      = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_R]      = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_G]      = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_oP1] = { // oP1
    .latt_type        = _oP1,
    .is_need_calculated = false,
    .n_highsym_points = 8,
    .n_paths          = 12,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_X,
      _hsp_Z,
      _hsp_U,
      _hsp_Y,
      _hsp_S,
      _hsp_T,
      _hsp_R, 
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_S},
      {_hsp_S,        _hsp_Y},
      {_hsp_Y,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_Z,        _hsp_U},
      {_hsp_U,        _hsp_R},
      {_hsp_R,        _hsp_T},
      {_hsp_T,        _hsp_Z},
      {_hsp_X,        _hsp_U},
      {_hsp_Y,        _hsp_T},
      {_hsp_S,        _hsp_R},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_X]      = {   0.5,     0,     0},
      [_hsp_Z]      = {     0,     0,   0.5},
      [_hsp_U]      = {   0.5,     0,   0.5},
      [_hsp_Y]      = {     0,   0.5,     0},
      [_hsp_S]      = {   0.5,   0.5,     0},
      [_hsp_T]      = {     0,   0.5,   0.5},
      [_hsp_R]      = {   0.5,   0.5,   0.5},
    },
  },
  [_oF1] = { // oF1
    .latt_type        = _oF1,
    .is_need_calculated = true,
    .n_highsym_points = 9,
    .n_paths          = 9,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_T,
      _hsp_Z,
      _hsp_Y,
      _hsp_SIGMA_0,
      _hsp_U_0,
      _hsp_A_0,
      _hsp_C_0,
      _hsp_L,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_Y},
      {_hsp_Y,        _hsp_T},
      {_hsp_T,        _hsp_Z},
      {_hsp_Z,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_SIGMA_0},
      {_hsp_U_0,      _hsp_T},
      {_hsp_Y,        _hsp_C_0},
      {_hsp_A_0,      _hsp_Z},
      {_hsp_GAMMA,    _hsp_L},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]    = {     0,     0,     0},
      [_hsp_T]        = {     1,   0.5,   0.5},
      [_hsp_Z]        = {   0.5,   0.5,     0},
      [_hsp_Y]        = {   0.5,     0,   0.5},
      [_hsp_SIGMA_0]  = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_U_0]      = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_A_0]      = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_C_0]      = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_L]        = {   0.5,   0.5,   0.5},
    },
  },
  [_oF2] = { // oF2 --> Need to be calculated
    .latt_type        = _oF2,
    .is_need_calculated = true,
    .n_highsym_points = 9,
    .n_paths          = 9,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_T,
      _hsp_Z,
      _hsp_Y,
      _hsp_LAMBDA_0,
      _hsp_Q_0,
      _hsp_G_0,
      _hsp_H_0,
      _hsp_L,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_T},
      {_hsp_T,        _hsp_Z},
      {_hsp_Z,        _hsp_Y},
      {_hsp_Y,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_LAMBDA_0},
      {_hsp_Q_0,      _hsp_Z},
      {_hsp_T,        _hsp_G_0},
      {_hsp_H_0,      _hsp_Y},
      {_hsp_GAMMA,    _hsp_L},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_T]      = {     0,   0.5,   0.5},
      [_hsp_Z]      = {   0.5,   0.5,     1},
      [_hsp_Y]      = {   0.5,     0,   0.5},
      [_hsp_LAMBDA_0] = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_Q_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_G_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_H_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_L]      = {   0.5,   0.5,   0.5},
    },
  },
  [_oF3] = { // oF3 --> Need to be calculated
    .latt_type        = _oF3,
    .is_need_calculated = true,
    .n_highsym_points = 11,
    .n_paths          = 10,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_T,
      _hsp_Z,
      _hsp_Y,
      _hsp_A_0,
      _hsp_C_0,
      _hsp_B_0,
      _hsp_D_0,
      _hsp_G_0,
      _hsp_H_0,
      _hsp_L,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_Y},
      {_hsp_Y,        _hsp_C_0},
      {_hsp_A_0,      _hsp_Z},
      {_hsp_Z,        _hsp_B_0},
      {_hsp_D_0,      _hsp_T},
      {_hsp_T,        _hsp_G_0},
      {_hsp_H_0,      _hsp_Y},
      {_hsp_T,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_GAMMA,    _hsp_L},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_T]      = {     0,   0.5,   0.5},
      [_hsp_Z]      = {   0.5,   0.5,     0},
      [_hsp_Y]      = {   0.5,     0,   0.5},
      [_hsp_A_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_C_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_B_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_D_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_G_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_H_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_L]      = {   0.5,   0.5,   0.5},
    },
  },
  [_oI1] = { // oI1 --> Need to be calculated
    .latt_type        = _oI1,
    .is_need_calculated = true,
    .n_highsym_points = 13,
    .n_paths          = 11,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_X,
      _hsp_S,
      _hsp_R,
      _hsp_T,
      _hsp_W,
      _hsp_SIGMA_0,
      _hsp_F_2,
      _hsp_Y_0,
      _hsp_U_0,
      _hsp_L_0,
      _hsp_M_0,
      _hsp_J_0,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_F_2},
      {_hsp_SIGMA_0,  _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Y_0},
      {_hsp_U_0,      _hsp_X},
      {_hsp_GAMMA,    _hsp_R},
      {_hsp_R,        _hsp_W},
      {_hsp_W,        _hsp_S},
      {_hsp_S,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_T},
      {_hsp_T,        _hsp_W},
    },
    .hsp_coordinates  = {
      [_hsp_GAMMA]  = {     0,     0,     0},
      [_hsp_X]      = {   0.5,   0.5,   -.5},
      [_hsp_S]      = {   0.5,     0,     0},
      [_hsp_R]      = {     0,   0.5,     0},
      [_hsp_T]      = {     0,     0,   0.5},
      [_hsp_W]      = {  0.25,  0.25,  0.25},
      [_hsp_SIGMA_0] = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_F_2]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_Y_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_U_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_L_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_M_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_J_0]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_oI2] = { // oI2 --> Need to be calculated
    .latt_type        = _oI2,
    .is_need_calculated = true,
    .n_highsym_points = 13,
    .n_paths          = 11,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_X,
     	_hsp_S,
     	_hsp_R,
     	_hsp_T,
     	_hsp_W,
     	_hsp_Y_0,
     	_hsp_U_2,
     	_hsp_LAMBDA_0,
     	_hsp_G_2,
     	_hsp_K,
     	_hsp_K_2,
     	_hsp_K_4,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_U_2},
      {_hsp_Y_0,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_LAMBDA_0},
      {_hsp_G_2,      _hsp_X},
      {_hsp_GAMMA,    _hsp_R},
      {_hsp_R,        _hsp_W},
      {_hsp_W,        _hsp_S},
      {_hsp_S,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_T},
      {_hsp_T,        _hsp_W},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_X]      = {   -.5,   0.5,   0.5},
     	[_hsp_S]      = {   0.5,     0,     0},
     	[_hsp_R]      = {     0,   0.5,     0},
     	[_hsp_T]      = {     0,     0,   0.5},
     	[_hsp_W]      = {  0.25,  0.25,  0.25},
     	[_hsp_Y_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_U_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_LAMBDA_0] = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_K]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_K_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_K_4]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_oI3] = { // oI3 --> Need to be calculated
    .latt_type        = _oI3,
    .is_need_calculated = true,
    .n_highsym_points = 13,
    .n_paths          = 11,
    .highsym_points   = {
      _hsp_GAMMA,
      _hsp_X,
      _hsp_S,
      _hsp_R,
      _hsp_T,
      _hsp_W,
      _hsp_SIGMA_0,
      _hsp_F_0,
      _hsp_LAMBDA_0,
      _hsp_G_0,
      _hsp_V_0,
      _hsp_H_0,
      _hsp_H_2,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_X,        _hsp_F_0},
      {_hsp_SIGMA_0,  _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_LAMBDA_0},
      {_hsp_G_0,      _hsp_X},
      {_hsp_GAMMA,    _hsp_R},
      {_hsp_R,        _hsp_W},
      {_hsp_W,        _hsp_S},
      {_hsp_S,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_T},
      {_hsp_T,        _hsp_W},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_X]      = {   -.5,   0.5,   0.5},
     	[_hsp_S]      = {   0.5,     0,     0},
     	[_hsp_R]      = {     0,   0.5,     0},
     	[_hsp_T]      = {     0,     0,   0.5},
     	[_hsp_W]      = {  0.25,  0.25,  0.25},
      [_hsp_SIGMA_0]  = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_F_0]      = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_LAMBDA_0] = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_G_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_V_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_H_0]    = {    -1,    -1,    -1}, // Need to be calculated
      [_hsp_H_2]    = {    -1,    -1,    -1}, // Need to be calculated
    }
  },
  [_oC1] = { // oC1 --> Need to be calculated
    .latt_type        = _oC1,
    .is_need_calculated = true,
    .n_highsym_points = 10,
    .n_paths          = 11,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Y,
     	_hsp_T,
     	_hsp_Z,
     	_hsp_S,
     	_hsp_R,
     	_hsp_SIGMA_0,
     	_hsp_C_0,
     	_hsp_A_0,
     	_hsp_E_0,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_Y},
      {_hsp_Y,        _hsp_C_0},
      {_hsp_SIGMA_0,  _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_Z,        _hsp_A_0},
      {_hsp_E_0,      _hsp_T},
      {_hsp_T,        _hsp_Y},
      {_hsp_GAMMA,    _hsp_S},
      {_hsp_S,        _hsp_R},
      {_hsp_R,        _hsp_Z},
      {_hsp_Z,        _hsp_T},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Y]      = {   -.5,   0.5,     0},
     	[_hsp_T]      = {   -.5,   0.5,   0.5},
     	[_hsp_Z]      = {     0,     0,   0.5},
     	[_hsp_S]      = {     0,   0.5,     0},
     	[_hsp_R]      = {     0,   0.5,   0.5},
     	[_hsp_SIGMA_0] = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_C_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_A_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_E_0]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_oC2] = { // oC2 --> Need to be calculated
    .latt_type        = _oC2,
    .is_need_calculated = true,
    .n_highsym_points = 15,
    .n_paths          = 11,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Y,
     	_hsp_T,
     	_hsp_T_2,
     	_hsp_Z,
     	_hsp_Z_2,
     	_hsp_S,
     	_hsp_R,
     	_hsp_R_2,
     	_hsp_DELTA_0,
     	_hsp_F_0,
     	_hsp_B_0,
     	_hsp_B_2,
     	_hsp_G_0,
     	_hsp_G_2,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_Y},
      {_hsp_Y,        _hsp_F_0},
      {_hsp_DELTA_0,  _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_Z,        _hsp_B_0},
      {_hsp_G_0,      _hsp_T},
      {_hsp_T,        _hsp_Y},
      {_hsp_GAMMA,    _hsp_S},
      {_hsp_S,        _hsp_R},
      {_hsp_R,        _hsp_Z},
      {_hsp_Z,        _hsp_T},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Y]      = {   0.5,   0.5,     0},
     	[_hsp_T]      = {   0.5,   0.5,   0.5},
     	[_hsp_T_2]    = {   0.5,   0.5,   -.5},
     	[_hsp_Z]      = {     0,     0,   0.5},
     	[_hsp_Z_2]    = {     0,     0,   -.5},
     	[_hsp_S]      = {     0,   0.5,     0},
     	[_hsp_R]      = {     0,     0,   0.5},
     	[_hsp_R_2]    = {     0,   0.5,   -.5},
     	[_hsp_DELTA_0] = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_F_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_B_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_B_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G_2]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_oA1] = { // oA1 --> Need to be calculated
    .latt_type        = _oA1,
    .is_need_calculated = true,
    .n_highsym_points = 10,
    .n_paths          = 11,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Y,
     	_hsp_T,
     	_hsp_Z,
     	_hsp_S,
     	_hsp_R,
     	_hsp_SIGMA_0,
     	_hsp_C_0,
     	_hsp_A_0,
     	_hsp_E_0,
    },
    .path = {
      {_hsp_GAMMA,    _hsp_Y},
      {_hsp_Y,        _hsp_C_0},
      {_hsp_SIGMA_0,  _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_Z,        _hsp_A_0},
      {_hsp_E_0,      _hsp_T},
      {_hsp_T,        _hsp_Y},
      {_hsp_GAMMA,    _hsp_S},
      {_hsp_S,        _hsp_R},
      {_hsp_R,        _hsp_Z},
      {_hsp_Z,        _hsp_T},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Y]      = {   -.5,   0.5,     0},
     	[_hsp_T]      = {   -.5,   0.5,   0.5},
     	[_hsp_Z]      = {     0,     0,   0.5},
     	[_hsp_S]      = {     0,   0.5,     0},
     	[_hsp_R]      = {     0,   0.5,   0.5},
     	[_hsp_SIGMA_0]  = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_C_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_A_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_E_0]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_oA2] = { // oA2 --> Need to be calculated
    .latt_type        = _oA2,
    .is_need_calculated = true,
    .n_highsym_points = 15,
    .n_paths          = 11,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Y,
     	_hsp_T,
     	_hsp_T_2,
     	_hsp_Z,
     	_hsp_Z_2,
     	_hsp_S,
     	_hsp_R,
     	_hsp_R_2,
     	_hsp_DELTA_0,
     	_hsp_F_0,
     	_hsp_B_0,
     	_hsp_B_2,
     	_hsp_G_0,
     	_hsp_G_2,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_Y},
      {_hsp_Y,        _hsp_F_0},
      {_hsp_DELTA_0,  _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_Z,        _hsp_B_0},
      {_hsp_G_0,      _hsp_T},
      {_hsp_T,        _hsp_Y},
      {_hsp_GAMMA,    _hsp_S},
      {_hsp_S,        _hsp_R},
      {_hsp_R,        _hsp_Z},
      {_hsp_Z,        _hsp_T},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Y]      = {   0.5,   0.5,     0},
     	[_hsp_T]      = {   0.5,   0.5,   0.5},
     	[_hsp_T_2]    = {   0.5,   0.5,   -.5},
     	[_hsp_Z]      = {     0,     0,   0.5},
     	[_hsp_Z_2]    = {     0,     0,   -.5},
     	[_hsp_S]      = {     0,   0.5,     0},
     	[_hsp_R]      = {     0,   0.5,   0.5},
     	[_hsp_R_2]    = {     0,   0.5,   -.5},
     	[_hsp_DELTA_0]  = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_F_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_B_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_B_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G_2]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_hP1] = { // hP1
    .latt_type        = _hP1,
    .is_need_calculated = false,
    .n_highsym_points = 7,
    .n_paths          = 10,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_A,
     	_hsp_K,
     	_hsp_H,
     	_hsp_H_2,
     	_hsp_M,
     	_hsp_L,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_M},
      {_hsp_M,        _hsp_K},
      {_hsp_K,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_A},
      {_hsp_A,        _hsp_L},
      {_hsp_L,        _hsp_H},
      {_hsp_H,        _hsp_A},
      {_hsp_L,        _hsp_M},
      {_hsp_H,        _hsp_K},
      {_hsp_K,        _hsp_H_2 },
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_A]      = {     0,     0,   0.5},
     	[_hsp_K]      = { 0.333, 0.333,     0},
     	[_hsp_H]      = { 0.333, 0.333,   0.5},
     	[_hsp_H_2]    = { 0.333, 0.333,   -.5},
     	[_hsp_M]      = {   0.5,     0,     0},
     	[_hsp_L]      = {   0.5,     0,   0.5},
    },
  },
  [_hP2] = { // hP2
    .latt_type        = _hP2,
    .is_need_calculated = false,
    .n_highsym_points = 7,
    .n_paths          = 9,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_A,
     	_hsp_K,
     	_hsp_H,
     	_hsp_H_2,
     	_hsp_M,
     	_hsp_L,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_M},
      {_hsp_M,        _hsp_K},
      {_hsp_K,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_A},
      {_hsp_A,        _hsp_L},
      {_hsp_L,        _hsp_H},
      {_hsp_H,        _hsp_A},
      {_hsp_L,        _hsp_M},
      {_hsp_H,        _hsp_K },
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_A]      = {     0,     0,   0.5},
     	[_hsp_K]      = { 0.333, 0.333,     0},
     	[_hsp_H]      = { 0.333, 0.333,   0.5},
     	[_hsp_H_2]    = { 0.333, 0.333,   -.5},
     	[_hsp_M]      = {   0.5,     0,     0},
     	[_hsp_L]      = {   0.5,     0,   0.5},
    },
  },
  [_hR1] = { // hR1 --> Need to be calculated
    .latt_type        = _hR1,
    .is_need_calculated = true,
    .n_highsym_points = 20,
    .n_paths          = 7,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_T,
     	_hsp_L,
     	_hsp_L_2,
     	_hsp_L_4,
     	_hsp_F,
     	_hsp_F_2,
     	_hsp_S_0,
     	_hsp_S_2,
     	_hsp_S_4,
     	_hsp_S_6,
     	_hsp_H_0,
     	_hsp_H_2,
     	_hsp_H_4,
     	_hsp_H_6,
     	_hsp_M_0,
     	_hsp_M_2,
     	_hsp_M_4,
     	_hsp_M_6,
     	_hsp_M_8,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_T},
      {_hsp_T,        _hsp_H_2},
      {_hsp_H_0,      _hsp_L},
      {_hsp_L,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_S_0},
      {_hsp_S_2,      _hsp_F},
      {_hsp_F,        _hsp_GAMMA },
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_T]      = {   0.5,   0.5,   0.5},
     	[_hsp_L]      = {   0.5,     0,     0},
     	[_hsp_L_2]    = {     0,   -.5,     0},
     	[_hsp_L_4]    = {     0,     0,   -.5},
     	[_hsp_F]      = {   0.5,     0,   0.5},
     	[_hsp_F_2]    = {   0.5,   0.5,     0},
     	[_hsp_S_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_S_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_S_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_S_6]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_6]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_M_0]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_M_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_M_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_M_6]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_M_8]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_hR2] = { // hR2 --> Need to be calculated
    .latt_type        = _hR2,
    .is_need_calculated = true,
    .n_highsym_points = 9,
    .n_paths          = 5,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_T,
     	_hsp_P_0,
     	_hsp_P_2,
     	_hsp_R_0,
     	_hsp_M,
     	_hsp_M_2,
     	_hsp_L,
     	_hsp_F,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_L},
      {_hsp_L,        _hsp_T},
      {_hsp_T,        _hsp_P_0},
      {_hsp_P_2,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_F},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_T]      = {   0.5,   -.5,   0.5},
     	[_hsp_P_0]    = {    -1,    -1,    -1},
     	[_hsp_P_2]    = {    -1,    -1,    -1},
     	[_hsp_R_0]    = {    -1,    -1,    -1},
     	[_hsp_M]      = {    -1,    -1,    -1},
     	[_hsp_M_2]    = {    -1,    -1,    -1},
     	[_hsp_L]      = {   0.5,     0,     0},
     	[_hsp_F]      = {   0.5,   -.5,     0},
    },
  },
  [_mP1] = { // mP1 --> Need to be calculated
    .latt_type        = _mP1,
    .is_need_calculated = true,
    .n_highsym_points = 18,
    .n_paths          = 10,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Z,
     	_hsp_B,
     	_hsp_B_2,
     	_hsp_Y,
     	_hsp_Y_2,
     	_hsp_C,
     	_hsp_C_2,
     	_hsp_D,
     	_hsp_D_2,
     	_hsp_A,
     	_hsp_E,
     	_hsp_H,
     	_hsp_H_2,
     	_hsp_H_4,
     	_hsp_M,
     	_hsp_M_2,
     	_hsp_M_4,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_Z,        _hsp_D},
      {_hsp_D,        _hsp_B},
      {_hsp_B,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_A},
      {_hsp_A,        _hsp_E},
      {_hsp_E,        _hsp_Z},
      {_hsp_Z,        _hsp_C_2},
      {_hsp_C_2,      _hsp_Y_2},
      {_hsp_Y_2,      _hsp_GAMMA},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Z]      = {     0,   0.5,     0},
     	[_hsp_B]      = {     0,     0,   0.5},
     	[_hsp_B_2]    = {     0,     0,   -.5},
     	[_hsp_Y]      = {   0.5,     0,     0},
     	[_hsp_Y_2]    = {   -.5,     0,     0},
     	[_hsp_C]      = {   0.5,   0.5,     0},
     	[_hsp_C_2]    = {   -.5,   0.5,     0},
     	[_hsp_D]      = {     0,   0.5,   0.5},
     	[_hsp_D_2]    = {     0,   0.5,   -.5},
     	[_hsp_A]      = {   -.5,     0,   0.5},
     	[_hsp_E]      = {   -.5,   0.5,   0.5},
     	[_hsp_H]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_M]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_M_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_M_4]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_mC1] = { // mC1 --> Need to be calculated
    .latt_type        = _mC1,
    .is_need_calculated = true,
    .n_highsym_points = 16,
    .n_paths          = 9,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Y_2,
     	_hsp_Y_4,
     	_hsp_A,
     	_hsp_M_2,
     	_hsp_V,
     	_hsp_V_2,
     	_hsp_L_2,
     	_hsp_C,
     	_hsp_C_2,
     	_hsp_C_4,
     	_hsp_D,
     	_hsp_D_2,
     	_hsp_E,
     	_hsp_E_2,
     	_hsp_E_4,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_C},
      {_hsp_C_2,      _hsp_Y_2},
      {_hsp_Y_2,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_M_2},
      {_hsp_M_2,      _hsp_D},
      {_hsp_D_2,      _hsp_A},
      {_hsp_A,        _hsp_GAMMA},
      {_hsp_L_2,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_V_2 },
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Y_2]    = {   -.5,   0.5,     0},
     	[_hsp_Y_4]    = {   0.5,   -.5,     0},
     	[_hsp_A]      = {     0,     0,   0.5},
     	[_hsp_M_2]    = {   -.5,   0.5,   0.5},
     	[_hsp_V]      = {   0.5,     0,     0},
     	[_hsp_V_2]    = {     0,   0.5,     0},
     	[_hsp_L_2]    = {     0,   0.5,   0.5},
     	[_hsp_C]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_C_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_C_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_D]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_D_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_E]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_E_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_E_4]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_mC2] = { // mC2 --> Need to be calculated
    .latt_type        = _mC2,
    .is_need_calculated = true,
    .n_highsym_points = 16,
    .n_paths          = 6,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Y,
     	_hsp_A,
     	_hsp_M,
     	_hsp_V_2,
     	_hsp_L_2,
     	_hsp_F,
     	_hsp_F_2,
     	_hsp_F_4,
     	_hsp_H,
     	_hsp_H_2,
     	_hsp_H_4,
     	_hsp_G,
     	_hsp_G_2,
     	_hsp_G_4,
     	_hsp_G_6,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_Y},
      {_hsp_Y,        _hsp_M},
      {_hsp_M,        _hsp_A},
      {_hsp_A,        _hsp_GAMMA},
      {_hsp_L_2,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_V_2},
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Y]      = {   0.5,   0.5,     0},
     	[_hsp_A]      = {     0,     0,   0.5},
     	[_hsp_M]      = {   0.5,   0.5,   0.5},
     	[_hsp_V_2]    = {     0,   0.5,     0},
     	[_hsp_L_2]    = {     0,   0.5,   0.5},
     	[_hsp_F]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_F_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_F_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_G_6]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_mC3] = { // mC3 --> Need to be calculated
    .latt_type        = _mC3,
    .is_need_calculated = true,
    .n_highsym_points = 19,
    .n_paths          = 7,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Y,
     	_hsp_A,
     	_hsp_M_2,
     	_hsp_V,
     	_hsp_V_2,
     	_hsp_L_2,
     	_hsp_I,
     	_hsp_I_2,
     	_hsp_K,
     	_hsp_K_2,
     	_hsp_K_4,
     	_hsp_H,
     	_hsp_H_2,
     	_hsp_H_4,
     	_hsp_N,
     	_hsp_N_2,
     	_hsp_N_4,
     	_hsp_N_6,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_A},
      {_hsp_A,        _hsp_I_2},
      {_hsp_I,        _hsp_M_2},
      {_hsp_M_2,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Y},
      {_hsp_L_2,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_V_2 },
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Y]      = {   0.5,   0.5,     0},
     	[_hsp_A]      = {     0,     0,   0.5},
     	[_hsp_M_2]    = {   -.5,   0.5,   0.5},
     	[_hsp_V]      = {   0.5,     0,     0},
     	[_hsp_V_2]    = {     0,   0.5,     0},
     	[_hsp_L_2]    = {     0,   0.5,   0.5},
     	[_hsp_I]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_I_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_K]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_K_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_K_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_H_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_N]      = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_N_2]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_N_4]    = {    -1,    -1,    -1}, // Need to be calculated
     	[_hsp_N_6]    = {    -1,    -1,    -1}, // Need to be calculated
    },
  },
  [_aP2] = { // aP2
    .latt_type        = _aP2,
    .is_need_calculated = false,
    .n_highsym_points = 8,
    .n_paths          = 7,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Z,
     	_hsp_Y,
     	_hsp_X,
     	_hsp_V,
     	_hsp_U,
     	_hsp_T,
     	_hsp_R,
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_Y,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_R,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_T},
      {_hsp_U,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_V },
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA] = {     0,     0,     0},
     	[_hsp_Z]     = {     0,     0,   0.5},
     	[_hsp_Y]     = {     0,   0.5,     0},
     	[_hsp_X]     = {   0.5,     0,     0},
     	[_hsp_V]     = {   0.5,   0.5,     0},
     	[_hsp_U]     = {   0.5,     0,   0.5},
     	[_hsp_T]     = {     0,   0.5,   0.5},
     	[_hsp_R]     = {   0.5,   0.5,   0.5},
    },
  },
  [_aP3] = { // aP3
    .latt_type        = _aP3,
    .is_need_calculated = false,
    .n_highsym_points = 9,
    .n_paths          = 7,
    .highsym_points   = {
     	_hsp_GAMMA,
     	_hsp_Z,
     	_hsp_Y,
     	_hsp_Y_2,
     	_hsp_X,
     	_hsp_V_2,
     	_hsp_U_2,
     	_hsp_T_2,
     	_hsp_R_2, 
    },
    .path             = {
      {_hsp_GAMMA,    _hsp_X},
      {_hsp_Y,        _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_Z},
      {_hsp_R_2,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_T_2},
      {_hsp_U_2,      _hsp_GAMMA},
      {_hsp_GAMMA,    _hsp_V_2 },
    },
    .hsp_coordinates  = {
     	[_hsp_GAMMA]  = {     0,     0,     0},
     	[_hsp_Z]      = {     0,     0,   0.5},
     	[_hsp_Y]      = {     0,   0.5,     0},
     	[_hsp_Y_2]    = {     0,   -.5,     0},
     	[_hsp_X]      = {   0.5,     0,     0},
     	[_hsp_V_2]    = {   0.5,   -.5,     0},
     	[_hsp_U_2]    = {   -.5,     0,   0.5},
     	[_hsp_T_2]    = {     0,   -.5,   0.5},
     	[_hsp_R_2]    = {   -.5,   -.5,   0.5},
    },
  },
};
