#ifndef _GAMMA_MULT_TABLE_H
#define _GAMMA_MULT_TABLE_H

namespace cvc {

#if 0
extern int gamma_mult_table[16][16];
extern double gamma_mult_sign[16][16];
extern double gamma_adjoint_sign[16];
extern double gamma_transposed_sign[16];
#endif  // of if 0

/***************************************************************************
 * the sequence is:
 * 0, 1, 2, 3, id, 5, 0_5, 1_5, 2_5, 3_5, 0_1, 0_2, 0_3, 1_2, 1_3, 2_3
 * 0  1  2  3   4  5    6    7    8    9   10   11   12   13   14   15
 ***************************************************************************/

inline void init_gamma_mult_table (const char *gamma_basis) {
if (strcmp(gamma_basis,"cvc")==0){

  g_gamma_mult_table[0][0] = 4;
  g_gamma_mult_sign[0][0]  = 1;
  g_gamma_mult_table[0][1] = 10;
  g_gamma_mult_sign[0][1]  = 1;
  g_gamma_mult_table[0][2] = 11;
  g_gamma_mult_sign[0][2]  = 1;
  g_gamma_mult_table[0][3] = 12;
  g_gamma_mult_sign[0][3]  = 1;
  g_gamma_mult_table[0][4] = 0;
  g_gamma_mult_sign[0][4]  = 1;
  g_gamma_mult_table[0][5] = 6;
  g_gamma_mult_sign[0][5]  = 1;
  g_gamma_mult_table[0][6] = 5;
  g_gamma_mult_sign[0][6]  = 1;
  g_gamma_mult_table[0][7] = 15;
  g_gamma_mult_sign[0][7]  = -1;
  g_gamma_mult_table[0][8] = 14;
  g_gamma_mult_sign[0][8]  = 1;
  g_gamma_mult_table[0][9] = 13;
  g_gamma_mult_sign[0][9]  = -1;
  g_gamma_mult_table[0][10] = 1;
  g_gamma_mult_sign[0][10]  = 1;
  g_gamma_mult_table[0][11] = 2;
  g_gamma_mult_sign[0][11]  = 1;
  g_gamma_mult_table[0][12] = 3;
  g_gamma_mult_sign[0][12]  = 1;
  g_gamma_mult_table[0][13] = 9;
  g_gamma_mult_sign[0][13]  = -1;
  g_gamma_mult_table[0][14] = 8;
  g_gamma_mult_sign[0][14]  = 1;
  g_gamma_mult_table[0][15] = 7;
  g_gamma_mult_sign[0][15]  = -1;
  g_gamma_mult_table[1][0] = 10;
  g_gamma_mult_sign[1][0]  = -1;
  g_gamma_mult_table[1][1] = 4;
  g_gamma_mult_sign[1][1]  = 1;
  g_gamma_mult_table[1][2] = 13;
  g_gamma_mult_sign[1][2]  = 1;
  g_gamma_mult_table[1][3] = 14;
  g_gamma_mult_sign[1][3]  = 1;
  g_gamma_mult_table[1][4] = 1;
  g_gamma_mult_sign[1][4]  = 1;
  g_gamma_mult_table[1][5] = 7;
  g_gamma_mult_sign[1][5]  = 1;
  g_gamma_mult_table[1][6] = 15;
  g_gamma_mult_sign[1][6]  = 1;
  g_gamma_mult_table[1][7] = 5;
  g_gamma_mult_sign[1][7]  = 1;
  g_gamma_mult_table[1][8] = 12;
  g_gamma_mult_sign[1][8]  = -1;
  g_gamma_mult_table[1][9] = 11;
  g_gamma_mult_sign[1][9]  = 1;
  g_gamma_mult_table[1][10] = 0;
  g_gamma_mult_sign[1][10]  = -1;
  g_gamma_mult_table[1][11] = 9;
  g_gamma_mult_sign[1][11]  = 1;
  g_gamma_mult_table[1][12] = 8;
  g_gamma_mult_sign[1][12]  = -1;
  g_gamma_mult_table[1][13] = 2;
  g_gamma_mult_sign[1][13]  = 1;
  g_gamma_mult_table[1][14] = 3;
  g_gamma_mult_sign[1][14]  = 1;
  g_gamma_mult_table[1][15] = 6;
  g_gamma_mult_sign[1][15]  = 1;
  g_gamma_mult_table[2][0] = 11;
  g_gamma_mult_sign[2][0]  = -1;
  g_gamma_mult_table[2][1] = 13;
  g_gamma_mult_sign[2][1]  = -1;
  g_gamma_mult_table[2][2] = 4;
  g_gamma_mult_sign[2][2]  = 1;
  g_gamma_mult_table[2][3] = 15;
  g_gamma_mult_sign[2][3]  = 1;
  g_gamma_mult_table[2][4] = 2;
  g_gamma_mult_sign[2][4]  = 1;
  g_gamma_mult_table[2][5] = 8;
  g_gamma_mult_sign[2][5]  = 1;
  g_gamma_mult_table[2][6] = 14;
  g_gamma_mult_sign[2][6]  = -1;
  g_gamma_mult_table[2][7] = 12;
  g_gamma_mult_sign[2][7]  = 1;
  g_gamma_mult_table[2][8] = 5;
  g_gamma_mult_sign[2][8]  = 1;
  g_gamma_mult_table[2][9] = 10;
  g_gamma_mult_sign[2][9]  = -1;
  g_gamma_mult_table[2][10] = 9;
  g_gamma_mult_sign[2][10]  = -1;
  g_gamma_mult_table[2][11] = 0;
  g_gamma_mult_sign[2][11]  = -1;
  g_gamma_mult_table[2][12] = 7;
  g_gamma_mult_sign[2][12]  = 1;
  g_gamma_mult_table[2][13] = 1;
  g_gamma_mult_sign[2][13]  = -1;
  g_gamma_mult_table[2][14] = 6;
  g_gamma_mult_sign[2][14]  = -1;
  g_gamma_mult_table[2][15] = 3;
  g_gamma_mult_sign[2][15]  = 1;
  g_gamma_mult_table[3][0] = 12;
  g_gamma_mult_sign[3][0]  = -1;
  g_gamma_mult_table[3][1] = 14;
  g_gamma_mult_sign[3][1]  = -1;
  g_gamma_mult_table[3][2] = 15;
  g_gamma_mult_sign[3][2]  = -1;
  g_gamma_mult_table[3][3] = 4;
  g_gamma_mult_sign[3][3]  = 1;
  g_gamma_mult_table[3][4] = 3;
  g_gamma_mult_sign[3][4]  = 1;
  g_gamma_mult_table[3][5] = 9;
  g_gamma_mult_sign[3][5]  = 1;
  g_gamma_mult_table[3][6] = 13;
  g_gamma_mult_sign[3][6]  = 1;
  g_gamma_mult_table[3][7] = 11;
  g_gamma_mult_sign[3][7]  = -1;
  g_gamma_mult_table[3][8] = 10;
  g_gamma_mult_sign[3][8]  = 1;
  g_gamma_mult_table[3][9] = 5;
  g_gamma_mult_sign[3][9]  = 1;
  g_gamma_mult_table[3][10] = 8;
  g_gamma_mult_sign[3][10]  = 1;
  g_gamma_mult_table[3][11] = 7;
  g_gamma_mult_sign[3][11]  = -1;
  g_gamma_mult_table[3][12] = 0;
  g_gamma_mult_sign[3][12]  = -1;
  g_gamma_mult_table[3][13] = 6;
  g_gamma_mult_sign[3][13]  = 1;
  g_gamma_mult_table[3][14] = 1;
  g_gamma_mult_sign[3][14]  = -1;
  g_gamma_mult_table[3][15] = 2;
  g_gamma_mult_sign[3][15]  = -1;
  g_gamma_mult_table[4][0] = 0;
  g_gamma_mult_sign[4][0]  = 1;
  g_gamma_mult_table[4][1] = 1;
  g_gamma_mult_sign[4][1]  = 1;
  g_gamma_mult_table[4][2] = 2;
  g_gamma_mult_sign[4][2]  = 1;
  g_gamma_mult_table[4][3] = 3;
  g_gamma_mult_sign[4][3]  = 1;
  g_gamma_mult_table[4][4] = 4;
  g_gamma_mult_sign[4][4]  = 1;
  g_gamma_mult_table[4][5] = 5;
  g_gamma_mult_sign[4][5]  = 1;
  g_gamma_mult_table[4][6] = 6;
  g_gamma_mult_sign[4][6]  = 1;
  g_gamma_mult_table[4][7] = 7;
  g_gamma_mult_sign[4][7]  = 1;
  g_gamma_mult_table[4][8] = 8;
  g_gamma_mult_sign[4][8]  = 1;
  g_gamma_mult_table[4][9] = 9;
  g_gamma_mult_sign[4][9]  = 1;
  g_gamma_mult_table[4][10] = 10;
  g_gamma_mult_sign[4][10]  = 1;
  g_gamma_mult_table[4][11] = 11;
  g_gamma_mult_sign[4][11]  = 1;
  g_gamma_mult_table[4][12] = 12;
  g_gamma_mult_sign[4][12]  = 1;
  g_gamma_mult_table[4][13] = 13;
  g_gamma_mult_sign[4][13]  = 1;
  g_gamma_mult_table[4][14] = 14;
  g_gamma_mult_sign[4][14]  = 1;
  g_gamma_mult_table[4][15] = 15;
  g_gamma_mult_sign[4][15]  = 1;
  g_gamma_mult_table[5][0] = 6;
  g_gamma_mult_sign[5][0]  = -1;
  g_gamma_mult_table[5][1] = 7;
  g_gamma_mult_sign[5][1]  = -1;
  g_gamma_mult_table[5][2] = 8;
  g_gamma_mult_sign[5][2]  = -1;
  g_gamma_mult_table[5][3] = 9;
  g_gamma_mult_sign[5][3]  = -1;
  g_gamma_mult_table[5][4] = 5;
  g_gamma_mult_sign[5][4]  = 1;
  g_gamma_mult_table[5][5] = 4;
  g_gamma_mult_sign[5][5]  = 1;
  g_gamma_mult_table[5][6] = 0;
  g_gamma_mult_sign[5][6]  = -1;
  g_gamma_mult_table[5][7] = 1;
  g_gamma_mult_sign[5][7]  = -1;
  g_gamma_mult_table[5][8] = 2;
  g_gamma_mult_sign[5][8]  = -1;
  g_gamma_mult_table[5][9] = 3;
  g_gamma_mult_sign[5][9]  = -1;
  g_gamma_mult_table[5][10] = 15;
  g_gamma_mult_sign[5][10]  = -1;
  g_gamma_mult_table[5][11] = 14;
  g_gamma_mult_sign[5][11]  = 1;
  g_gamma_mult_table[5][12] = 13;
  g_gamma_mult_sign[5][12]  = -1;
  g_gamma_mult_table[5][13] = 12;
  g_gamma_mult_sign[5][13]  = -1;
  g_gamma_mult_table[5][14] = 11;
  g_gamma_mult_sign[5][14]  = 1;
  g_gamma_mult_table[5][15] = 10;
  g_gamma_mult_sign[5][15]  = -1;
  g_gamma_mult_table[6][0] = 5;
  g_gamma_mult_sign[6][0]  = -1;
  g_gamma_mult_table[6][1] = 15;
  g_gamma_mult_sign[6][1]  = 1;
  g_gamma_mult_table[6][2] = 14;
  g_gamma_mult_sign[6][2]  = -1;
  g_gamma_mult_table[6][3] = 13;
  g_gamma_mult_sign[6][3]  = 1;
  g_gamma_mult_table[6][4] = 6;
  g_gamma_mult_sign[6][4]  = 1;
  g_gamma_mult_table[6][5] = 0;
  g_gamma_mult_sign[6][5]  = 1;
  g_gamma_mult_table[6][6] = 4;
  g_gamma_mult_sign[6][6]  = -1;
  g_gamma_mult_table[6][7] = 10;
  g_gamma_mult_sign[6][7]  = -1;
  g_gamma_mult_table[6][8] = 11;
  g_gamma_mult_sign[6][8]  = -1;
  g_gamma_mult_table[6][9] = 12;
  g_gamma_mult_sign[6][9]  = -1;
  g_gamma_mult_table[6][10] = 7;
  g_gamma_mult_sign[6][10]  = 1;
  g_gamma_mult_table[6][11] = 8;
  g_gamma_mult_sign[6][11]  = 1;
  g_gamma_mult_table[6][12] = 9;
  g_gamma_mult_sign[6][12]  = 1;
  g_gamma_mult_table[6][13] = 3;
  g_gamma_mult_sign[6][13]  = -1;
  g_gamma_mult_table[6][14] = 2;
  g_gamma_mult_sign[6][14]  = 1;
  g_gamma_mult_table[6][15] = 1;
  g_gamma_mult_sign[6][15]  = -1;
  g_gamma_mult_table[7][0] = 15;
  g_gamma_mult_sign[7][0]  = -1;
  g_gamma_mult_table[7][1] = 5;
  g_gamma_mult_sign[7][1]  = -1;
  g_gamma_mult_table[7][2] = 12;
  g_gamma_mult_sign[7][2]  = 1;
  g_gamma_mult_table[7][3] = 11;
  g_gamma_mult_sign[7][3]  = -1;
  g_gamma_mult_table[7][4] = 7;
  g_gamma_mult_sign[7][4]  = 1;
  g_gamma_mult_table[7][5] = 1;
  g_gamma_mult_sign[7][5]  = 1;
  g_gamma_mult_table[7][6] = 10;
  g_gamma_mult_sign[7][6]  = 1;
  g_gamma_mult_table[7][7] = 4;
  g_gamma_mult_sign[7][7]  = -1;
  g_gamma_mult_table[7][8] = 13;
  g_gamma_mult_sign[7][8]  = -1;
  g_gamma_mult_table[7][9] = 14;
  g_gamma_mult_sign[7][9]  = -1;
  g_gamma_mult_table[7][10] = 6;
  g_gamma_mult_sign[7][10]  = -1;
  g_gamma_mult_table[7][11] = 3;
  g_gamma_mult_sign[7][11]  = 1;
  g_gamma_mult_table[7][12] = 2;
  g_gamma_mult_sign[7][12]  = -1;
  g_gamma_mult_table[7][13] = 8;
  g_gamma_mult_sign[7][13]  = 1;
  g_gamma_mult_table[7][14] = 9;
  g_gamma_mult_sign[7][14]  = 1;
  g_gamma_mult_table[7][15] = 0;
  g_gamma_mult_sign[7][15]  = 1;
  g_gamma_mult_table[8][0] = 14;
  g_gamma_mult_sign[8][0]  = 1;
  g_gamma_mult_table[8][1] = 12;
  g_gamma_mult_sign[8][1]  = -1;
  g_gamma_mult_table[8][2] = 5;
  g_gamma_mult_sign[8][2]  = -1;
  g_gamma_mult_table[8][3] = 10;
  g_gamma_mult_sign[8][3]  = 1;
  g_gamma_mult_table[8][4] = 8;
  g_gamma_mult_sign[8][4]  = 1;
  g_gamma_mult_table[8][5] = 2;
  g_gamma_mult_sign[8][5]  = 1;
  g_gamma_mult_table[8][6] = 11;
  g_gamma_mult_sign[8][6]  = 1;
  g_gamma_mult_table[8][7] = 13;
  g_gamma_mult_sign[8][7]  = 1;
  g_gamma_mult_table[8][8] = 4;
  g_gamma_mult_sign[8][8]  = -1;
  g_gamma_mult_table[8][9] = 15;
  g_gamma_mult_sign[8][9]  = -1;
  g_gamma_mult_table[8][10] = 3;
  g_gamma_mult_sign[8][10]  = -1;
  g_gamma_mult_table[8][11] = 6;
  g_gamma_mult_sign[8][11]  = -1;
  g_gamma_mult_table[8][12] = 1;
  g_gamma_mult_sign[8][12]  = 1;
  g_gamma_mult_table[8][13] = 7;
  g_gamma_mult_sign[8][13]  = -1;
  g_gamma_mult_table[8][14] = 0;
  g_gamma_mult_sign[8][14]  = -1;
  g_gamma_mult_table[8][15] = 9;
  g_gamma_mult_sign[8][15]  = 1;
  g_gamma_mult_table[9][0] = 13;
  g_gamma_mult_sign[9][0]  = -1;
  g_gamma_mult_table[9][1] = 11;
  g_gamma_mult_sign[9][1]  = 1;
  g_gamma_mult_table[9][2] = 10;
  g_gamma_mult_sign[9][2]  = -1;
  g_gamma_mult_table[9][3] = 5;
  g_gamma_mult_sign[9][3]  = -1;
  g_gamma_mult_table[9][4] = 9;
  g_gamma_mult_sign[9][4]  = 1;
  g_gamma_mult_table[9][5] = 3;
  g_gamma_mult_sign[9][5]  = 1;
  g_gamma_mult_table[9][6] = 12;
  g_gamma_mult_sign[9][6]  = 1;
  g_gamma_mult_table[9][7] = 14;
  g_gamma_mult_sign[9][7]  = 1;
  g_gamma_mult_table[9][8] = 15;
  g_gamma_mult_sign[9][8]  = 1;
  g_gamma_mult_table[9][9] = 4;
  g_gamma_mult_sign[9][9]  = -1;
  g_gamma_mult_table[9][10] = 2;
  g_gamma_mult_sign[9][10]  = 1;
  g_gamma_mult_table[9][11] = 1;
  g_gamma_mult_sign[9][11]  = -1;
  g_gamma_mult_table[9][12] = 6;
  g_gamma_mult_sign[9][12]  = -1;
  g_gamma_mult_table[9][13] = 0;
  g_gamma_mult_sign[9][13]  = 1;
  g_gamma_mult_table[9][14] = 7;
  g_gamma_mult_sign[9][14]  = -1;
  g_gamma_mult_table[9][15] = 8;
  g_gamma_mult_sign[9][15]  = -1;
  g_gamma_mult_table[10][0] = 1;
  g_gamma_mult_sign[10][0]  = -1;
  g_gamma_mult_table[10][1] = 0;
  g_gamma_mult_sign[10][1]  = 1;
  g_gamma_mult_table[10][2] = 9;
  g_gamma_mult_sign[10][2]  = -1;
  g_gamma_mult_table[10][3] = 8;
  g_gamma_mult_sign[10][3]  = 1;
  g_gamma_mult_table[10][4] = 10;
  g_gamma_mult_sign[10][4]  = 1;
  g_gamma_mult_table[10][5] = 15;
  g_gamma_mult_sign[10][5]  = -1;
  g_gamma_mult_table[10][6] = 7;
  g_gamma_mult_sign[10][6]  = -1;
  g_gamma_mult_table[10][7] = 6;
  g_gamma_mult_sign[10][7]  = 1;
  g_gamma_mult_table[10][8] = 3;
  g_gamma_mult_sign[10][8]  = -1;
  g_gamma_mult_table[10][9] = 2;
  g_gamma_mult_sign[10][9]  = 1;
  g_gamma_mult_table[10][10] = 4;
  g_gamma_mult_sign[10][10]  = -1;
  g_gamma_mult_table[10][11] = 13;
  g_gamma_mult_sign[10][11]  = -1;
  g_gamma_mult_table[10][12] = 14;
  g_gamma_mult_sign[10][12]  = -1;
  g_gamma_mult_table[10][13] = 11;
  g_gamma_mult_sign[10][13]  = 1;
  g_gamma_mult_table[10][14] = 12;
  g_gamma_mult_sign[10][14]  = 1;
  g_gamma_mult_table[10][15] = 5;
  g_gamma_mult_sign[10][15]  = 1;
  g_gamma_mult_table[11][0] = 2;
  g_gamma_mult_sign[11][0]  = -1;
  g_gamma_mult_table[11][1] = 9;
  g_gamma_mult_sign[11][1]  = 1;
  g_gamma_mult_table[11][2] = 0;
  g_gamma_mult_sign[11][2]  = 1;
  g_gamma_mult_table[11][3] = 7;
  g_gamma_mult_sign[11][3]  = -1;
  g_gamma_mult_table[11][4] = 11;
  g_gamma_mult_sign[11][4]  = 1;
  g_gamma_mult_table[11][5] = 14;
  g_gamma_mult_sign[11][5]  = 1;
  g_gamma_mult_table[11][6] = 8;
  g_gamma_mult_sign[11][6]  = -1;
  g_gamma_mult_table[11][7] = 3;
  g_gamma_mult_sign[11][7]  = 1;
  g_gamma_mult_table[11][8] = 6;
  g_gamma_mult_sign[11][8]  = 1;
  g_gamma_mult_table[11][9] = 1;
  g_gamma_mult_sign[11][9]  = -1;
  g_gamma_mult_table[11][10] = 13;
  g_gamma_mult_sign[11][10]  = 1;
  g_gamma_mult_table[11][11] = 4;
  g_gamma_mult_sign[11][11]  = -1;
  g_gamma_mult_table[11][12] = 15;
  g_gamma_mult_sign[11][12]  = -1;
  g_gamma_mult_table[11][13] = 10;
  g_gamma_mult_sign[11][13]  = -1;
  g_gamma_mult_table[11][14] = 5;
  g_gamma_mult_sign[11][14]  = -1;
  g_gamma_mult_table[11][15] = 12;
  g_gamma_mult_sign[11][15]  = 1;
  g_gamma_mult_table[12][0] = 3;
  g_gamma_mult_sign[12][0]  = -1;
  g_gamma_mult_table[12][1] = 8;
  g_gamma_mult_sign[12][1]  = -1;
  g_gamma_mult_table[12][2] = 7;
  g_gamma_mult_sign[12][2]  = 1;
  g_gamma_mult_table[12][3] = 0;
  g_gamma_mult_sign[12][3]  = 1;
  g_gamma_mult_table[12][4] = 12;
  g_gamma_mult_sign[12][4]  = 1;
  g_gamma_mult_table[12][5] = 13;
  g_gamma_mult_sign[12][5]  = -1;
  g_gamma_mult_table[12][6] = 9;
  g_gamma_mult_sign[12][6]  = -1;
  g_gamma_mult_table[12][7] = 2;
  g_gamma_mult_sign[12][7]  = -1;
  g_gamma_mult_table[12][8] = 1;
  g_gamma_mult_sign[12][8]  = 1;
  g_gamma_mult_table[12][9] = 6;
  g_gamma_mult_sign[12][9]  = 1;
  g_gamma_mult_table[12][10] = 14;
  g_gamma_mult_sign[12][10]  = 1;
  g_gamma_mult_table[12][11] = 15;
  g_gamma_mult_sign[12][11]  = 1;
  g_gamma_mult_table[12][12] = 4;
  g_gamma_mult_sign[12][12]  = -1;
  g_gamma_mult_table[12][13] = 5;
  g_gamma_mult_sign[12][13]  = 1;
  g_gamma_mult_table[12][14] = 10;
  g_gamma_mult_sign[12][14]  = -1;
  g_gamma_mult_table[12][15] = 11;
  g_gamma_mult_sign[12][15]  = -1;
  g_gamma_mult_table[13][0] = 9;
  g_gamma_mult_sign[13][0]  = -1;
  g_gamma_mult_table[13][1] = 2;
  g_gamma_mult_sign[13][1]  = -1;
  g_gamma_mult_table[13][2] = 1;
  g_gamma_mult_sign[13][2]  = 1;
  g_gamma_mult_table[13][3] = 6;
  g_gamma_mult_sign[13][3]  = 1;
  g_gamma_mult_table[13][4] = 13;
  g_gamma_mult_sign[13][4]  = 1;
  g_gamma_mult_table[13][5] = 12;
  g_gamma_mult_sign[13][5]  = -1;
  g_gamma_mult_table[13][6] = 3;
  g_gamma_mult_sign[13][6]  = -1;
  g_gamma_mult_table[13][7] = 8;
  g_gamma_mult_sign[13][7]  = -1;
  g_gamma_mult_table[13][8] = 7;
  g_gamma_mult_sign[13][8]  = 1;
  g_gamma_mult_table[13][9] = 0;
  g_gamma_mult_sign[13][9]  = 1;
  g_gamma_mult_table[13][10] = 11;
  g_gamma_mult_sign[13][10]  = -1;
  g_gamma_mult_table[13][11] = 10;
  g_gamma_mult_sign[13][11]  = 1;
  g_gamma_mult_table[13][12] = 5;
  g_gamma_mult_sign[13][12]  = 1;
  g_gamma_mult_table[13][13] = 4;
  g_gamma_mult_sign[13][13]  = -1;
  g_gamma_mult_table[13][14] = 15;
  g_gamma_mult_sign[13][14]  = -1;
  g_gamma_mult_table[13][15] = 14;
  g_gamma_mult_sign[13][15]  = 1;
  g_gamma_mult_table[14][0] = 8;
  g_gamma_mult_sign[14][0]  = 1;
  g_gamma_mult_table[14][1] = 3;
  g_gamma_mult_sign[14][1]  = -1;
  g_gamma_mult_table[14][2] = 6;
  g_gamma_mult_sign[14][2]  = -1;
  g_gamma_mult_table[14][3] = 1;
  g_gamma_mult_sign[14][3]  = 1;
  g_gamma_mult_table[14][4] = 14;
  g_gamma_mult_sign[14][4]  = 1;
  g_gamma_mult_table[14][5] = 11;
  g_gamma_mult_sign[14][5]  = 1;
  g_gamma_mult_table[14][6] = 2;
  g_gamma_mult_sign[14][6]  = 1;
  g_gamma_mult_table[14][7] = 9;
  g_gamma_mult_sign[14][7]  = -1;
  g_gamma_mult_table[14][8] = 0;
  g_gamma_mult_sign[14][8]  = -1;
  g_gamma_mult_table[14][9] = 7;
  g_gamma_mult_sign[14][9]  = 1;
  g_gamma_mult_table[14][10] = 12;
  g_gamma_mult_sign[14][10]  = -1;
  g_gamma_mult_table[14][11] = 5;
  g_gamma_mult_sign[14][11]  = -1;
  g_gamma_mult_table[14][12] = 10;
  g_gamma_mult_sign[14][12]  = 1;
  g_gamma_mult_table[14][13] = 15;
  g_gamma_mult_sign[14][13]  = 1;
  g_gamma_mult_table[14][14] = 4;
  g_gamma_mult_sign[14][14]  = -1;
  g_gamma_mult_table[14][15] = 13;
  g_gamma_mult_sign[14][15]  = -1;
  g_gamma_mult_table[15][0] = 7;
  g_gamma_mult_sign[15][0]  = -1;
  g_gamma_mult_table[15][1] = 6;
  g_gamma_mult_sign[15][1]  = 1;
  g_gamma_mult_table[15][2] = 3;
  g_gamma_mult_sign[15][2]  = -1;
  g_gamma_mult_table[15][3] = 2;
  g_gamma_mult_sign[15][3]  = 1;
  g_gamma_mult_table[15][4] = 15;
  g_gamma_mult_sign[15][4]  = 1;
  g_gamma_mult_table[15][5] = 10;
  g_gamma_mult_sign[15][5]  = -1;
  g_gamma_mult_table[15][6] = 1;
  g_gamma_mult_sign[15][6]  = -1;
  g_gamma_mult_table[15][7] = 0;
  g_gamma_mult_sign[15][7]  = 1;
  g_gamma_mult_table[15][8] = 9;
  g_gamma_mult_sign[15][8]  = -1;
  g_gamma_mult_table[15][9] = 8;
  g_gamma_mult_sign[15][9]  = 1;
  g_gamma_mult_table[15][10] = 5;
  g_gamma_mult_sign[15][10]  = 1;
  g_gamma_mult_table[15][11] = 12;
  g_gamma_mult_sign[15][11]  = -1;
  g_gamma_mult_table[15][12] = 11;
  g_gamma_mult_sign[15][12]  = 1;
  g_gamma_mult_table[15][13] = 14;
  g_gamma_mult_sign[15][13]  = -1;
  g_gamma_mult_table[15][14] = 13;
  g_gamma_mult_sign[15][14]  = 1;
  g_gamma_mult_table[15][15] = 4;
  g_gamma_mult_sign[15][15]  = -1;


  g_gamma_adjoint_sign[0] = 1;
  g_gamma_adjoint_sign[1] = 1;
  g_gamma_adjoint_sign[2] = 1;
  g_gamma_adjoint_sign[3] = 1;
  g_gamma_adjoint_sign[4] = 1;
  g_gamma_adjoint_sign[5] = 1;
  g_gamma_adjoint_sign[6] = -1;
  g_gamma_adjoint_sign[7] = -1;
  g_gamma_adjoint_sign[8] = -1;
  g_gamma_adjoint_sign[9] = -1;
  g_gamma_adjoint_sign[10] = -1;
  g_gamma_adjoint_sign[11] = -1;
  g_gamma_adjoint_sign[12] = -1;
  g_gamma_adjoint_sign[13] = -1;
  g_gamma_adjoint_sign[14] = -1;
  g_gamma_adjoint_sign[15] = -1;



  g_gamma_transposed_sign[0] = 1;
  g_gamma_transposed_sign[1] = -1;
  g_gamma_transposed_sign[2] = 1;
  g_gamma_transposed_sign[3] = -1;
  g_gamma_transposed_sign[4] = 1;
  g_gamma_transposed_sign[5] = 1;
  g_gamma_transposed_sign[6] = -1;
  g_gamma_transposed_sign[7] = 1;
  g_gamma_transposed_sign[8] = -1;
  g_gamma_transposed_sign[9] = 1;
  g_gamma_transposed_sign[10] = 1;
  g_gamma_transposed_sign[11] = -1;
  g_gamma_transposed_sign[12] = 1;
  g_gamma_transposed_sign[13] = 1;
  g_gamma_transposed_sign[14] = -1;
  g_gamma_transposed_sign[15] = 1;
} else if (strcmp(gamma_basis,"plegma")==0){
  g_gamma_mult_table[0][0] = 4;
  g_gamma_mult_sign[0][0]  = 1;
  g_gamma_mult_table[0][1] = 10;
  g_gamma_mult_sign[0][1]  = 1;
  g_gamma_mult_table[0][2] = 11;
  g_gamma_mult_sign[0][2]  = 1;
  g_gamma_mult_table[0][3] = 12;
  g_gamma_mult_sign[0][3]  = 1;
  g_gamma_mult_table[0][4] = 0;
  g_gamma_mult_sign[0][4]  = 1;
  g_gamma_mult_table[0][5] = 6;
  g_gamma_mult_sign[0][5]  = 1;
  g_gamma_mult_table[0][6] = 5;
  g_gamma_mult_sign[0][6]  = 1;
  g_gamma_mult_table[0][7] = 15;
  g_gamma_mult_sign[0][7]  = 1;
  g_gamma_mult_table[0][8] = 14;
  g_gamma_mult_sign[0][8]  = -1;
  g_gamma_mult_table[0][9] = 13;
  g_gamma_mult_sign[0][9]  = 1;
  g_gamma_mult_table[0][10] = 1;
  g_gamma_mult_sign[0][10]  = 1;
  g_gamma_mult_table[0][11] = 2;
  g_gamma_mult_sign[0][11]  = 1;
  g_gamma_mult_table[0][12] = 3;
  g_gamma_mult_sign[0][12]  = 1;
  g_gamma_mult_table[0][13] = 9;
  g_gamma_mult_sign[0][13]  = 1;
  g_gamma_mult_table[0][14] = 8;
  g_gamma_mult_sign[0][14]  = -1;
  g_gamma_mult_table[0][15] = 7;
  g_gamma_mult_sign[0][15]  = 1;

  g_gamma_mult_table[1][0] = 10;
  g_gamma_mult_sign[1][0]  = -1;
  g_gamma_mult_table[1][1] = 4;
  g_gamma_mult_sign[1][1]  = 1;
  g_gamma_mult_table[1][2] = 13;
  g_gamma_mult_sign[1][2]  = 1;
  g_gamma_mult_table[1][3] = 14;
  g_gamma_mult_sign[1][3]  = 1;
  g_gamma_mult_table[1][4] = 1;
  g_gamma_mult_sign[1][4]  = 1;
  g_gamma_mult_table[1][5] = 7;
  g_gamma_mult_sign[1][5]  = 1;
  g_gamma_mult_table[1][6] = 15;
  g_gamma_mult_sign[1][6]  = -1;
  g_gamma_mult_table[1][7] = 5;
  g_gamma_mult_sign[1][7]  = 1;
  g_gamma_mult_table[1][8] = 12;
  g_gamma_mult_sign[1][8]  = 1;
  g_gamma_mult_table[1][9] = 11;
  g_gamma_mult_sign[1][9]  = -1;
  g_gamma_mult_table[1][10] = 0;
  g_gamma_mult_sign[1][10]  = -1;
  g_gamma_mult_table[1][11] = 9;
  g_gamma_mult_sign[1][11]  = -1;
  g_gamma_mult_table[1][12] = 8;
  g_gamma_mult_sign[1][12]  = 1;
  g_gamma_mult_table[1][13] = 2;
  g_gamma_mult_sign[1][13]  = 1;
  g_gamma_mult_table[1][14] = 3;
  g_gamma_mult_sign[1][14]  = 1;
  g_gamma_mult_table[1][15] = 6;
  g_gamma_mult_sign[1][15]  = -1;


  g_gamma_mult_table[2][0] = 11;
  g_gamma_mult_sign[2][0]  = -1;
  g_gamma_mult_table[2][1] = 13;
  g_gamma_mult_sign[2][1]  = -1;
  g_gamma_mult_table[2][2] = 4;
  g_gamma_mult_sign[2][2]  = 1;
  g_gamma_mult_table[2][3] = 15;
  g_gamma_mult_sign[2][3]  = 1;
  g_gamma_mult_table[2][4] = 2;
  g_gamma_mult_sign[2][4]  = 1;
  g_gamma_mult_table[2][5] = 8;
  g_gamma_mult_sign[2][5]  = 1;
  g_gamma_mult_table[2][6] = 14;
  g_gamma_mult_sign[2][6]  = 1;
  g_gamma_mult_table[2][7] = 12;
  g_gamma_mult_sign[2][7]  = -1;
  g_gamma_mult_table[2][8] = 5;
  g_gamma_mult_sign[2][8]  = 1;
  g_gamma_mult_table[2][9] = 10;
  g_gamma_mult_sign[2][9]  = 1;
  g_gamma_mult_table[2][10] = 9;
  g_gamma_mult_sign[2][10]  = 1;
  g_gamma_mult_table[2][11] = 0;
  g_gamma_mult_sign[2][11]  = -1;
  g_gamma_mult_table[2][12] = 7;
  g_gamma_mult_sign[2][12]  = -1;
  g_gamma_mult_table[2][13] = 1;
  g_gamma_mult_sign[2][13]  = -1;
  g_gamma_mult_table[2][14] = 6;
  g_gamma_mult_sign[2][14]  = 1;
  g_gamma_mult_table[2][15] = 3;
  g_gamma_mult_sign[2][15]  = 1;

  g_gamma_mult_table[3][0] = 12;
  g_gamma_mult_sign[3][0]  = -1;
  g_gamma_mult_table[3][1] = 14;
  g_gamma_mult_sign[3][1]  = -1;
  g_gamma_mult_table[3][2] = 15;
  g_gamma_mult_sign[3][2]  = -1;
  g_gamma_mult_table[3][3] = 4;
  g_gamma_mult_sign[3][3]  = 1;
  g_gamma_mult_table[3][4] = 3;
  g_gamma_mult_sign[3][4]  = 1;
  g_gamma_mult_table[3][5] = 9;
  g_gamma_mult_sign[3][5]  = 1;
  g_gamma_mult_table[3][6] = 13;
  g_gamma_mult_sign[3][6]  = 1;
  g_gamma_mult_table[3][7] = 11;
  g_gamma_mult_sign[3][7]  = 1;
  g_gamma_mult_table[3][8] = 10;
  g_gamma_mult_sign[3][8]  = -1;
  g_gamma_mult_table[3][9] = 5;
  g_gamma_mult_sign[3][9]  = 1;
  g_gamma_mult_table[3][10] = 8;
  g_gamma_mult_sign[3][10]  = -1;
  g_gamma_mult_table[3][11] = 7;
  g_gamma_mult_sign[3][11]  = 1;
  g_gamma_mult_table[3][12] = 0;
  g_gamma_mult_sign[3][12]  = -1;
  g_gamma_mult_table[3][13] = 6;
  g_gamma_mult_sign[3][13]  = 1;
  g_gamma_mult_table[3][14] = 1;
  g_gamma_mult_sign[3][14]  = -1;
  g_gamma_mult_table[3][15] = 2;
  g_gamma_mult_sign[3][15]  = -1;

  g_gamma_mult_table[4][0] = 0;
  g_gamma_mult_sign[4][0]  = 1;
  g_gamma_mult_table[4][1] = 1;
  g_gamma_mult_sign[4][1]  = 1;
  g_gamma_mult_table[4][2] = 2;
  g_gamma_mult_sign[4][2]  = 1;
  g_gamma_mult_table[4][3] = 3;
  g_gamma_mult_sign[4][3]  = 1;
  g_gamma_mult_table[4][4] = 4;
  g_gamma_mult_sign[4][4]  = 1;
  g_gamma_mult_table[4][5] = 5;
  g_gamma_mult_sign[4][5]  = 1;
  g_gamma_mult_table[4][6] = 6;
  g_gamma_mult_sign[4][6]  = 1;
  g_gamma_mult_table[4][7] = 7;
  g_gamma_mult_sign[4][7]  = 1;
  g_gamma_mult_table[4][8] = 8;
  g_gamma_mult_sign[4][8]  = 1;
  g_gamma_mult_table[4][9] = 9;
  g_gamma_mult_sign[4][9]  = 1;
  g_gamma_mult_table[4][10] = 10;
  g_gamma_mult_sign[4][10]  = 1;
  g_gamma_mult_table[4][11] = 11;
  g_gamma_mult_sign[4][11]  = 1;
  g_gamma_mult_table[4][12] = 12;
  g_gamma_mult_sign[4][12]  = 1;
  g_gamma_mult_table[4][13] = 13;
  g_gamma_mult_sign[4][13]  = 1;
  g_gamma_mult_table[4][14] = 14;
  g_gamma_mult_sign[4][14]  = 1;
  g_gamma_mult_table[4][15] = 15;
  g_gamma_mult_sign[4][15]  = 1;

  g_gamma_mult_table[5][0] = 6;
  g_gamma_mult_sign[5][0]  = -1;
  g_gamma_mult_table[5][1] = 7;
  g_gamma_mult_sign[5][1]  = -1;
  g_gamma_mult_table[5][2] = 8;
  g_gamma_mult_sign[5][2]  = -1;
  g_gamma_mult_table[5][3] = 9;
  g_gamma_mult_sign[5][3]  = -1;
  g_gamma_mult_table[5][4] = 5;
  g_gamma_mult_sign[5][4]  = 1;
  g_gamma_mult_table[5][5] = 4;
  g_gamma_mult_sign[5][5]  = 1;
  g_gamma_mult_table[5][6] = 0;
  g_gamma_mult_sign[5][6]  = -1;
  g_gamma_mult_table[5][7] = 1;
  g_gamma_mult_sign[5][7]  = -1;
  g_gamma_mult_table[5][8] = 2;
  g_gamma_mult_sign[5][8]  = -1;
  g_gamma_mult_table[5][9] = 3;
  g_gamma_mult_sign[5][9]  = -1;
  g_gamma_mult_table[5][10] = 15;
  g_gamma_mult_sign[5][10]  = 1;
  g_gamma_mult_table[5][11] = 14;
  g_gamma_mult_sign[5][11]  = -1;
  g_gamma_mult_table[5][12] = 13;
  g_gamma_mult_sign[5][12]  =  1;
  g_gamma_mult_table[5][13] = 12;
  g_gamma_mult_sign[5][13]  = -1;
  g_gamma_mult_table[5][14] = 11;
  g_gamma_mult_sign[5][14]  = -1;
  g_gamma_mult_table[5][15] = 10;
  g_gamma_mult_sign[5][15]  = +1;

  g_gamma_mult_table[6][0] = 5;
  g_gamma_mult_sign[6][0]  = -1;
  g_gamma_mult_table[6][1] = 15;
  g_gamma_mult_sign[6][1]  = -1;
  g_gamma_mult_table[6][2] = 14;
  g_gamma_mult_sign[6][2]  = 1;
  g_gamma_mult_table[6][3] = 13;
  g_gamma_mult_sign[6][3]  = -1;
  g_gamma_mult_table[6][4] = 6;
  g_gamma_mult_sign[6][4]  = 1;
  g_gamma_mult_table[6][5] = 0;
  g_gamma_mult_sign[6][5]  = 1;
  g_gamma_mult_table[6][6] = 4;
  g_gamma_mult_sign[6][6]  = -1;
  g_gamma_mult_table[6][7] = 10;
  g_gamma_mult_sign[6][7]  = -1;
  g_gamma_mult_table[6][8] = 11;
  g_gamma_mult_sign[6][8]  = -1;
  g_gamma_mult_table[6][9] = 12;
  g_gamma_mult_sign[6][9]  = -1;
  g_gamma_mult_table[6][10] = 7;
  g_gamma_mult_sign[6][10]  = 1;
  g_gamma_mult_table[6][11] = 8;
  g_gamma_mult_sign[6][11]  = 1;
  g_gamma_mult_table[6][12] = 9;
  g_gamma_mult_sign[6][12]  = 1;
  g_gamma_mult_table[6][13] = 3;
  g_gamma_mult_sign[6][13]  = 1;
  g_gamma_mult_table[6][14] = 2;
  g_gamma_mult_sign[6][14]  = -1;
  g_gamma_mult_table[6][15] = 1;
  g_gamma_mult_sign[6][15]  = 1;

  g_gamma_mult_table[7][0] = 15;
  g_gamma_mult_sign[7][0]  = 1;
  g_gamma_mult_table[7][1] = 5;
  g_gamma_mult_sign[7][1]  = -1;
  g_gamma_mult_table[7][2] = 12;
  g_gamma_mult_sign[7][2]  = -1;
  g_gamma_mult_table[7][3] = 11;
  g_gamma_mult_sign[7][3]  = 1;
  g_gamma_mult_table[7][4] = 7;
  g_gamma_mult_sign[7][4]  = 1;
  g_gamma_mult_table[7][5] = 1;
  g_gamma_mult_sign[7][5]  = 1;
  g_gamma_mult_table[7][6] = 10;
  g_gamma_mult_sign[7][6]  = 1;
  g_gamma_mult_table[7][7] = 4;
  g_gamma_mult_sign[7][7]  = -1;
  g_gamma_mult_table[7][8] = 13;
  g_gamma_mult_sign[7][8]  = -1;
  g_gamma_mult_table[7][9] = 14;
  g_gamma_mult_sign[7][9]  = -1;
  g_gamma_mult_table[7][10] = 6;
  g_gamma_mult_sign[7][10]  = -1;
  g_gamma_mult_table[7][11] = 3;
  g_gamma_mult_sign[7][11]  = -1;
  g_gamma_mult_table[7][12] = 2;
  g_gamma_mult_sign[7][12]  = 1;
  g_gamma_mult_table[7][13] = 8;
  g_gamma_mult_sign[7][13]  = 1;
  g_gamma_mult_table[7][14] = 9;
  g_gamma_mult_sign[7][14]  = 1;
  g_gamma_mult_table[7][15] = 0;
  g_gamma_mult_sign[7][15]  = -1;

  g_gamma_mult_table[8][0] = 14;
  g_gamma_mult_sign[8][0]  = -1;
  g_gamma_mult_table[8][1] = 12;
  g_gamma_mult_sign[8][1]  = 1;
  g_gamma_mult_table[8][2] = 5;
  g_gamma_mult_sign[8][2]  = -1;
  g_gamma_mult_table[8][3] = 10;
  g_gamma_mult_sign[8][3]  = -1;
  g_gamma_mult_table[8][4] = 8;
  g_gamma_mult_sign[8][4]  = 1;
  g_gamma_mult_table[8][5] = 2;
  g_gamma_mult_sign[8][5]  = 1;
  g_gamma_mult_table[8][6] = 11;
  g_gamma_mult_sign[8][6]  = 1;
  g_gamma_mult_table[8][7] = 13;
  g_gamma_mult_sign[8][7]  = 1;
  g_gamma_mult_table[8][8] = 4;
  g_gamma_mult_sign[8][8]  = -1;
  g_gamma_mult_table[8][9] = 15;
  g_gamma_mult_sign[8][9]  = -1;
  g_gamma_mult_table[8][10] = 3;
  g_gamma_mult_sign[8][10]  = 1;
  g_gamma_mult_table[8][11] = 6;
  g_gamma_mult_sign[8][11]  = -1;
  g_gamma_mult_table[8][12] = 1;
  g_gamma_mult_sign[8][12]  = -1;
  g_gamma_mult_table[8][13] = 7;
  g_gamma_mult_sign[8][13]  = -1;
  g_gamma_mult_table[8][14] = 0;
  g_gamma_mult_sign[8][14]  = 1;
  g_gamma_mult_table[8][15] = 9;
  g_gamma_mult_sign[8][15]  = 1;

  g_gamma_mult_table[9][0] = 13;
  g_gamma_mult_sign[9][0]  = 1;
  g_gamma_mult_table[9][1] = 11;
  g_gamma_mult_sign[9][1]  = -1;
  g_gamma_mult_table[9][2] = 10;
  g_gamma_mult_sign[9][2]  = 1;
  g_gamma_mult_table[9][3] = 5;
  g_gamma_mult_sign[9][3]  = -1;
  g_gamma_mult_table[9][4] = 9;
  g_gamma_mult_sign[9][4]  = 1;
  g_gamma_mult_table[9][5] = 3;
  g_gamma_mult_sign[9][5]  = 1;
  g_gamma_mult_table[9][6] = 12;
  g_gamma_mult_sign[9][6]  = 1;
  g_gamma_mult_table[9][7] = 14;
  g_gamma_mult_sign[9][7]  = 1;
  g_gamma_mult_table[9][8] = 15;
  g_gamma_mult_sign[9][8]  = 1;
  g_gamma_mult_table[9][9] = 4;
  g_gamma_mult_sign[9][9]  = -1;
  g_gamma_mult_table[9][10] = 2;
  g_gamma_mult_sign[9][10]  = -1;
  g_gamma_mult_table[9][11] = 1;
  g_gamma_mult_sign[9][11]  = 1;
  g_gamma_mult_table[9][12] = 6;
  g_gamma_mult_sign[9][12]  = -1;
  g_gamma_mult_table[9][13] = 0;
  g_gamma_mult_sign[9][13]  = -1;
  g_gamma_mult_table[9][14] = 7;
  g_gamma_mult_sign[9][14]  = -1;
  g_gamma_mult_table[9][15] = 8;
  g_gamma_mult_sign[9][15]  = -1;

  g_gamma_mult_table[10][0] = 1;
  g_gamma_mult_sign[10][0]  = -1;
  g_gamma_mult_table[10][1] = 0;
  g_gamma_mult_sign[10][1]  = 1;
  g_gamma_mult_table[10][2] = 9;
  g_gamma_mult_sign[10][2]  = 1;
  g_gamma_mult_table[10][3] = 8;
  g_gamma_mult_sign[10][3]  = -1;
  g_gamma_mult_table[10][4] = 10;
  g_gamma_mult_sign[10][4]  = 1;
  g_gamma_mult_table[10][5] = 15;
  g_gamma_mult_sign[10][5]  = 1;
  g_gamma_mult_table[10][6] = 7;
  g_gamma_mult_sign[10][6]  = -1;
  g_gamma_mult_table[10][7] = 6;
  g_gamma_mult_sign[10][7]  = 1;
  g_gamma_mult_table[10][8] = 3;
  g_gamma_mult_sign[10][8]  = -1;
  g_gamma_mult_table[10][9] = 2;
  g_gamma_mult_sign[10][9]  = -1;
  g_gamma_mult_table[10][10] = 4;
  g_gamma_mult_sign[10][10]  = -1;
  g_gamma_mult_table[10][11] = 13;
  g_gamma_mult_sign[10][11]  = -1;
  g_gamma_mult_table[10][12] = 14;
  g_gamma_mult_sign[10][12]  = -1;
  g_gamma_mult_table[10][13] = 11;
  g_gamma_mult_sign[10][13]  = 1;
  g_gamma_mult_table[10][14] = 12;
  g_gamma_mult_sign[10][14]  = 1;
  g_gamma_mult_table[10][15] = 5;
  g_gamma_mult_sign[10][15]  = -1;

  g_gamma_mult_table[11][0] = 2;
  g_gamma_mult_sign[11][0]  = -1;
  g_gamma_mult_table[11][1] = 9;
  g_gamma_mult_sign[11][1]  = -1;
  g_gamma_mult_table[11][2] = 0;
  g_gamma_mult_sign[11][2]  = 1;
  g_gamma_mult_table[11][3] = 7;
  g_gamma_mult_sign[11][3]  = 1;
  g_gamma_mult_table[11][4] = 11;
  g_gamma_mult_sign[11][4]  = 1;
  g_gamma_mult_table[11][5] = 14;
  g_gamma_mult_sign[11][5]  = -1;
  g_gamma_mult_table[11][6] = 8;
  g_gamma_mult_sign[11][6]  = -1;
  g_gamma_mult_table[11][7] = 3;
  g_gamma_mult_sign[11][7]  = -1;
  g_gamma_mult_table[11][8] = 6;
  g_gamma_mult_sign[11][8]  = 1;
  g_gamma_mult_table[11][9] = 1;
  g_gamma_mult_sign[11][9]  = 1;
  g_gamma_mult_table[11][10] = 13;
  g_gamma_mult_sign[11][10]  = 1;
  g_gamma_mult_table[11][11] = 4;
  g_gamma_mult_sign[11][11]  = -1;
  g_gamma_mult_table[11][12] = 15;
  g_gamma_mult_sign[11][12]  = -1;
  g_gamma_mult_table[11][13] = 10;
  g_gamma_mult_sign[11][13]  = -1;
  g_gamma_mult_table[11][14] = 5;
  g_gamma_mult_sign[11][14]  = 1;
  g_gamma_mult_table[11][15] = 12;
  g_gamma_mult_sign[11][15]  = 1;

  g_gamma_mult_table[12][0] = 3;
  g_gamma_mult_sign[12][0]  = -1;
  g_gamma_mult_table[12][1] = 8;
  g_gamma_mult_sign[12][1]  = 1;
  g_gamma_mult_table[12][2] = 7;
  g_gamma_mult_sign[12][2]  = -1;
  g_gamma_mult_table[12][3] = 0;
  g_gamma_mult_sign[12][3]  = 1;
  g_gamma_mult_table[12][4] = 12;
  g_gamma_mult_sign[12][4]  = 1;
  g_gamma_mult_table[12][5] = 13;
  g_gamma_mult_sign[12][5]  = 1;
  g_gamma_mult_table[12][6] = 9;
  g_gamma_mult_sign[12][6]  = -1;
  g_gamma_mult_table[12][7] = 2;
  g_gamma_mult_sign[12][7]  = 1;
  g_gamma_mult_table[12][8] = 1;
  g_gamma_mult_sign[12][8]  = -1;
  g_gamma_mult_table[12][9] = 6;
  g_gamma_mult_sign[12][9]  = 1;
  g_gamma_mult_table[12][10] = 14;
  g_gamma_mult_sign[12][10]  = 1;
  g_gamma_mult_table[12][11] = 15;
  g_gamma_mult_sign[12][11]  = 1;
  g_gamma_mult_table[12][12] = 4;
  g_gamma_mult_sign[12][12]  = -1;
  g_gamma_mult_table[12][13] = 5;
  g_gamma_mult_sign[12][13]  = -1;
  g_gamma_mult_table[12][14] = 10;
  g_gamma_mult_sign[12][14]  = -1;
  g_gamma_mult_table[12][15] = 11;
  g_gamma_mult_sign[12][15]  = -1;

  g_gamma_mult_table[13][0] = 9;
  g_gamma_mult_sign[13][0]  = 1;
  g_gamma_mult_table[13][1] = 2;
  g_gamma_mult_sign[13][1]  = -1;
  g_gamma_mult_table[13][2] = 1;
  g_gamma_mult_sign[13][2]  = 1;
  g_gamma_mult_table[13][3] = 6;
  g_gamma_mult_sign[13][3]  = -1;
  g_gamma_mult_table[13][4] = 13;
  g_gamma_mult_sign[13][4]  = 1;
  g_gamma_mult_table[13][5] = 12;
  g_gamma_mult_sign[13][5]  = -1;
  g_gamma_mult_table[13][6] = 3;
  g_gamma_mult_sign[13][6]  = 1;
  g_gamma_mult_table[13][7] = 8;
  g_gamma_mult_sign[13][7]  = -1;
  g_gamma_mult_table[13][8] = 7;
  g_gamma_mult_sign[13][8]  = 1;
  g_gamma_mult_table[13][9] = 0;
  g_gamma_mult_sign[13][9]  =-1;
  g_gamma_mult_table[13][10] = 11;
  g_gamma_mult_sign[13][10]  = -1;
  g_gamma_mult_table[13][11] = 10;
  g_gamma_mult_sign[13][11]  = 1;
  g_gamma_mult_table[13][12] = 5;
  g_gamma_mult_sign[13][12]  = -1;
  g_gamma_mult_table[13][13] = 4;
  g_gamma_mult_sign[13][13]  = -1;
  g_gamma_mult_table[13][14] = 15;
  g_gamma_mult_sign[13][14]  = -1;
  g_gamma_mult_table[13][15] = 14;
  g_gamma_mult_sign[13][15]  = 1;

  g_gamma_mult_table[14][0] = 8;
  g_gamma_mult_sign[14][0]  = -1;
  g_gamma_mult_table[14][1] = 3;
  g_gamma_mult_sign[14][1]  = -1;
  g_gamma_mult_table[14][2] = 6;
  g_gamma_mult_sign[14][2]  = 1;
  g_gamma_mult_table[14][3] = 1;
  g_gamma_mult_sign[14][3]  = 1;
  g_gamma_mult_table[14][4] = 14;
  g_gamma_mult_sign[14][4]  = 1;
  g_gamma_mult_table[14][5] = 11;
  g_gamma_mult_sign[14][5]  = -1;
  g_gamma_mult_table[14][6] = 2;
  g_gamma_mult_sign[14][6]  = -1;
  g_gamma_mult_table[14][7] = 9;
  g_gamma_mult_sign[14][7]  = -1;
  g_gamma_mult_table[14][8] = 0;
  g_gamma_mult_sign[14][8]  = 1;
  g_gamma_mult_table[14][9] = 7;
  g_gamma_mult_sign[14][9]  = 1;
  g_gamma_mult_table[14][10] = 12;
  g_gamma_mult_sign[14][10]  = -1;
  g_gamma_mult_table[14][11] = 5;
  g_gamma_mult_sign[14][11]  = 1;
  g_gamma_mult_table[14][12] = 10;
  g_gamma_mult_sign[14][12]  = 1;
  g_gamma_mult_table[14][13] = 15;
  g_gamma_mult_sign[14][13]  = 1;
  g_gamma_mult_table[14][14] = 4;
  g_gamma_mult_sign[14][14]  = -1;
  g_gamma_mult_table[14][15] = 13;
  g_gamma_mult_sign[14][15]  = -1;

  g_gamma_mult_table[15][0] = 7;
  g_gamma_mult_sign[15][0]  = 1;
  g_gamma_mult_table[15][1] = 6;
  g_gamma_mult_sign[15][1]  = -1;
  g_gamma_mult_table[15][2] = 3;
  g_gamma_mult_sign[15][2]  = -1;
  g_gamma_mult_table[15][3] = 2;
  g_gamma_mult_sign[15][3]  = 1;
  g_gamma_mult_table[15][4] = 15;
  g_gamma_mult_sign[15][4]  = 1;
  g_gamma_mult_table[15][5] = 10;
  g_gamma_mult_sign[15][5]  = 1;
  g_gamma_mult_table[15][6] = 1;
  g_gamma_mult_sign[15][6]  = 1;
  g_gamma_mult_table[15][7] = 0;
  g_gamma_mult_sign[15][7]  = -1;
  g_gamma_mult_table[15][8] = 9;
  g_gamma_mult_sign[15][8]  = -1;
  g_gamma_mult_table[15][9] = 8;
  g_gamma_mult_sign[15][9]  = 1;
  g_gamma_mult_table[15][10] = 5;
  g_gamma_mult_sign[15][10]  = -1;
  g_gamma_mult_table[15][11] = 12;
  g_gamma_mult_sign[15][11]  = -1;
  g_gamma_mult_table[15][12] = 11;
  g_gamma_mult_sign[15][12]  = 1;
  g_gamma_mult_table[15][13] = 14;
  g_gamma_mult_sign[15][13]  = -1;
  g_gamma_mult_table[15][14] = 13;
  g_gamma_mult_sign[15][14]  = 1;
  g_gamma_mult_table[15][15] = 4;
  g_gamma_mult_sign[15][15]  = -1;

  g_gamma_adjoint_sign[0] = 1;
  g_gamma_adjoint_sign[1] = 1;
  g_gamma_adjoint_sign[2] = 1;
  g_gamma_adjoint_sign[3] = 1;
  g_gamma_adjoint_sign[4] = 1;
  g_gamma_adjoint_sign[5] = 1;
  g_gamma_adjoint_sign[6] = -1;
  g_gamma_adjoint_sign[7] = -1;
  g_gamma_adjoint_sign[8] = -1;
  g_gamma_adjoint_sign[9] = -1;
  g_gamma_adjoint_sign[10] = -1;
  g_gamma_adjoint_sign[11] = -1;
  g_gamma_adjoint_sign[12] = -1;
  g_gamma_adjoint_sign[13] = -1;
  g_gamma_adjoint_sign[14] = -1;
  g_gamma_adjoint_sign[15] = -1;

  g_gamma_transposed_sign[0] = 1;
  g_gamma_transposed_sign[1] = -1;
  g_gamma_transposed_sign[2] = 1;
  g_gamma_transposed_sign[3] = -1;
  g_gamma_transposed_sign[4] = 1;
  g_gamma_transposed_sign[5] = 1;
  g_gamma_transposed_sign[6] = -1;
  g_gamma_transposed_sign[7] = 1;
  g_gamma_transposed_sign[8] = -1;
  g_gamma_transposed_sign[9] = 1;
  g_gamma_transposed_sign[10] = 1;
  g_gamma_transposed_sign[11] = -1;
  g_gamma_transposed_sign[12] = 1;
  g_gamma_transposed_sign[13] = 1;
  g_gamma_transposed_sign[14] = -1;
  g_gamma_transposed_sign[15] = 1;


}
}  // end of init_gamma_mult_table
}  // end of namespace cvc
#endif
