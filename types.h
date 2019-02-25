#ifndef _TYPES_H
#define _TYPES_H

#include <map>
#include <vector>
#include <string>
#include <list>

namespace cvc {

typedef double * spinor_vector_type;
typedef double * fermion_vector_type;
typedef double ** fermion_propagator_type;
typedef double ** spinor_propagator_type;

typedef struct mom_t {
  int x;
  int y;
  int z;
} mom_t;

typedef struct {
  int gi;
  int gf;
  int pi[3];
  int pf[3];
} m_m_2pt_type;

typedef struct {
  int gi1;
  int gi2;
  int gf1;
  int gf2;
  int pi1[3];
  int pi2[3];
  int pf1[3];
  int pf2[3];
} mxb_mxb_2pt_type;

typedef struct {
  double *****v;
  double *ptr;
  int np;
  int ng;
  int nt;
  int nv;
} gsp_type;

typedef std::map<std::string, std::vector< mom_t > > mom_lists_t;

typedef struct H5Correlator {
  // exists and does nothing to make this default-constructible, but
  // in practice the other constructor should be used
  H5Correlator() {}
  
  // note that 'storage', being a vector, will be zero-initialised
  H5Correlator(const std::list<std::string> & path_list_in,
               const size_t nelem) :
    path_list(path_list_in), storage(nelem) {}

  std::list<std::string> path_list;
  std::vector<double> storage;
} H5Correlator;

typedef struct deriv_t {
  int dim;
  int dir; 
} deriv_t;

}  /* end of namespace cvc */
#endif
