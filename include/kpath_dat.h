#pragma once

#include "error.h"
#include "spg_wrap.h"


#include "spglib/spg_database.h"
#include "spglib/spacegroup.h"    // enum Centering

/*
 * Explaination of the extended Bravais lattice symbol:
 *                ``xY#_''
 *  Where, `x' -> lattice system, including 
 *      `c' -> Cubic
 *      `t' -> Tetragonal
 *      `o' -> Orthorhombic
 *      `h' -> Hexagonal & Rhombohedral
 *      `m' -> Monoclinic
 *      `a' -> Triclinic
 *    `Y' -> Certering type, including
 *      `P' -> Primitive
 *      `F' -> Face centered
 *      `I' -> Body centered
 *      `S' -> Side-face centered, also known as `C' and `A'
 *      `R' -> Triple hexagonal
 *    `#' -> Index to discriminate the lattcies in same xY type
 *
 *  In spglib, SpglibDataset->
 */
typedef enum {
  _null = 0,
  _cP1,  /* 195 <= N <= 206 */
  _cP2,  /* 207 <= N <= 230 */
  _cF1,  /* 195 <= N <= 206 */
  _cF2,  /* 207 <= N <= 230 */
  _cI1,
  _tP1,
  _tI1,
  _tI2,
  _oP1,
  _oF1,
  _oF2,
  _oF3,
  _oI1,
  _oI2,
  _oI3,
  _oC1,
  _oC2,
  _oA1,
  _oA2,
  _hP1,
  _hP2,
  _hR1,
  _hR2,
  _mP1,
  _mC1,
  _mC2,
  _mC3,
  _aP2,
  _aP3,
} LatticeType;


LatticeType kpt_get_lattice_type(const Cellp* cell,
                                 const double prec,
                                 const int    hall_number);

ERROR_CODE kpt_print_lattice_type(const LatticeType latt_type,
                                        char        res_str[4]);




/*
 * k path related stuff
 */

typedef enum {
         _hsp_A = 0,  _hsp_A_0,
         _hsp_B,      _hsp_B_0,      _hsp_B_2,
         _hsp_C,      _hsp_C_0,      _hsp_C_2,      _hsp_C_4,
         _hsp_D,      _hsp_D_0,      _hsp_D_2,      
                  _hsp_DELTA_0,                          
         _hsp_E,      _hsp_E_0,      _hsp_E_2,      _hsp_E_4,
         _hsp_F,      _hsp_F_0,      _hsp_F_2,      _hsp_F_4,
         _hsp_G,      _hsp_G_0,      _hsp_G_2,      _hsp_G_4,      _hsp_G_6,
     _hsp_GAMMA,                                    
         _hsp_H,      _hsp_H_0,      _hsp_H_2,      _hsp_H_4,      _hsp_H_6,
         _hsp_I,                     _hsp_I_2,                          
                      _hsp_J_0,                                    
         _hsp_K,                     _hsp_K_2,      _hsp_K_4,                
         _hsp_L,      _hsp_L_0,      _hsp_L_2,      _hsp_L_4,      
  _hsp_LAMBDA_0,                                    
         _hsp_M,      _hsp_M_0,      _hsp_M_2,      _hsp_M_4,      _hsp_M_6,     _hsp_M_8,
         _hsp_N,      _hsp_N_2,      _hsp_N_4,      _hsp_N_6,      
         _hsp_P,      _hsp_P_0,      _hsp_P_2,                
                      _hsp_Q_0,                                    
         _hsp_R,      _hsp_R_0,      _hsp_R_2,                
         _hsp_S,      _hsp_S_0,      _hsp_S_2,      _hsp_S_4,      _hsp_S_6,
                  _hsp_SIGMA_0,                          
         _hsp_T,                     _hsp_T_2,                
         _hsp_U,      _hsp_U_0,      _hsp_U_2,      
         _hsp_V,      _hsp_V_0,      _hsp_V_2,      
         _hsp_W,                     _hsp_W_2,                
         _hsp_X,                     _hsp_X_1,                
         _hsp_Y,      _hsp_Y_0,      _hsp_Y_2,      _hsp_Y_4,
         _hsp_Z,      _hsp_Z_0,      _hsp_Z_2,
} HighSymmetryPoint;


typedef struct {
  LatticeType       latt_type;
  bool              is_need_calculated;
  int               n_highsym_points;
  int               n_paths;
  HighSymmetryPoint highsym_points[20];
  HighSymmetryPoint path[12][2];        // Maximum path nodes is 12
  double            hsp_coordinates[85][3];
} KPath;


KPath kpt_get_kpath(const LatticeType latt_type);
