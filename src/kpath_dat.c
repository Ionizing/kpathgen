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

