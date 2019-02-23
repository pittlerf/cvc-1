#pragma once

#include "cvc_complex.h"
#include "types.h"
#include "DependencyGraph.hpp"

#include <map>
#include <vector>
#include <string>
#include <cstring>
#include <stdexcept>

namespace cvc {

typedef struct ts_stoch_src_meta_t
{
  mom_t p;
  int gamma;
  int src_ts;
  
  ts_stoch_src_meta_t() {}
  
  ts_stoch_src_meta_t(const mom_t & p_in,
                      const int gamma_in,
                      const int src_ts_in) :
   p(p_in), gamma(gamma_in), src_ts(src_ts_in) {}

  std::string key(void)
  {
    return key(p, gamma, src_ts);
  }

  static std::string key(const mom_t & p_in,
                         const int gamma_in,
                         const int src_ts_in)
  {
    char key[100];
    snprintf(key, 100, "g%d/px%dpy%dpz%d/t%d", 
             gamma_in, 
             p_in.x, p_in.y, p_in.z, 
             src_ts_in);
    return std::string(key);
  }
} ts_stoch_src_meta_t;

static inline std::string check_mom_prop(const std::string prop_in)
{
  if( !(prop_in == "fwd" || prop_in == "bwd" || prop_in == "neither") ){
    throw std::invalid_argument("propagator with momentum has to be either 'fwd', 'bwd' or 'neither'\n");
  }
  return prop_in;
}

typedef struct stoch_prop_meta_t
{
  mom_t p;
  int gamma;
  int src_ts;
  std::string flav;
  std::string solver_driver;
  int solver_id;

  // this intentionally exists but does nothing
  stoch_prop_meta_t() {}
  
  stoch_prop_meta_t(const mom_t & p_in,
                    const int gamma_in,
                    const int src_ts_in,
                    const std::string & flav_in) :
    p(p_in), gamma(gamma_in), flav(flav_in), src_ts(src_ts_in),
    solver_driver("undefined"), solver_id(-1) {}

  stoch_prop_meta_t(const mom_t & p_in,
                    const int gamma_in,
                    const int src_ts_in,
                    const std::string & flav_in,
                    const std::string & solver_driver_in,
                    const int solver_id_in) :
    p(p_in), gamma(gamma_in), flav(flav_in), src_ts(src_ts_in),
    solver_driver(solver_driver_in), solver_id(solver_id_in) {}

  std::string key(void) const 
  {
    return key(p, gamma, src_ts, flav);
  }

  static std::string key(const mom_t & p_in,
                         const int gamma_in,
                         const int src_ts_in,
                         const std::string & flav_in)
  {
    char temp[100];
    snprintf(temp, 100, 
             "S%s_g%d_px%dpy%dpz%d_t%d",
             flav_in.c_str(), gamma_in, p_in.x, p_in.y, p_in.z, src_ts_in);
    return std::string(temp);
  }

} stoch_prop_meta_t;

typedef struct seq_stoch_prop_meta_t
{
  // this intentionally exists but does nothing
  seq_stoch_prop_meta_t() {}

  seq_stoch_prop_meta_t(const mom_t & seq_p_in,
                        const int seq_gamma_in,
                        const int seq_src_ts_in,
                        const int src_ts_in,
                        const std::string seq_flav_in,
                        const mom_t & src_p_in,
                        const int src_gamma_in,
                        const std::string & src_flav_in)
    : src_prop(src_p_in, src_gamma_in, src_ts_in, src_flav_in),
      p(seq_p_in),
      seq_src_ts(seq_src_ts_in),
      gamma(seq_gamma_in),
      flav(seq_flav_in),
      src_ts(src_ts_in) {}

  std::string key(void) const
  {
    return key(p, gamma, seq_src_ts, flav);
  }
  
  std::string key(const mom_t & p_in,
                  const int gamma_in,
                  const int seq_src_ts_in,
                  const std::string & flav_in) const
  {
    char temp[100];
    snprintf(temp, 100,
             "S%s_g%d_px%dpy%dpz%d::ts%d::",
             flav_in.c_str(), gamma_in, p_in.x, p_in.y, p_in.z, seq_src_ts_in);
    return( (std::string(temp)+src_prop.key()) );
  }

  mom_t p;
  int gamma;
  int src_ts;
  int seq_src_ts;
  std::string flav;

  stoch_prop_meta_t src_prop;
} seq_stoch_prop_meta_t;

//typedef struct shifted_prop_meta_t
//{
//  std::vector<shift_t> shifts;
//  std::string prop_key;
//} shifted_prop_meta_t;

/**
 * @brief Meta-description of a meson two-point function
 *
 * Describes the connected part of the two-point function
 *
 * \sum_{x_i, x_f}  e^{-i p_i x_i} e^{i p_f x_f}
 * \bar{chi}_{f_1}(x_f) gf chi_{f_2}(x_f) \bar{chi}_{f_2}(x_i) gi f_1(x_i)
 *
 * in the twisted basis. This is computed as
 *
 * Tr[ g5 (S_{f_1})^dag g5 gf S_{f_2} gi ]
 * = -[D^{-1}_{f_1^dag} gb eta]^dag g5 gf [D^{-1}_{f_2} g_i eta]  (1)
 *
 * where (f_1)^dag is bprop_flav (so if f_1 = d, bprop_flav = u)
 *        f_2      is fprop_flav
 *
 * and 'gb' in this case is just 'g5', but see comment below. Note that
 * the g5 multiplying gf is implicitly included in the later contraction.
 *
 * The momentum at source can either be included in the backward propagator
 * (in which case the phase is complex conjugated) or the forward propagator,
 * such that one can specify either "fwd" or "bwd" for src_mom_prop
 *
 * Similarly, in the writing of (1), it can be that it is more convenient to change
 * the gamma structure in the backward propagator, so we can specify 'gb' freely with the caveat
 * that what happens under complex conjugation must be considered externally,
 * unlike for the momentum.
 *
 */

typedef struct oet_meson_twopt_meta_t {
  oet_meson_twopt_meta_t(const std::string fprop_flav_in,
                   const std::string bprop_flav_in,
                   const std::string src_mom_prop_in,
                   const int gi_in,
                   const int gf_in,
                   const int gb_in,
                   const ::cvc::complex normalisation_in)
  {
    gi = gi_in;
    gf = gf_in;
    gb = gb_in;
    fprop_flav = fprop_flav_in;
    bprop_flav = bprop_flav_in;
    src_mom_prop = check_mom_prop(src_mom_prop_in);
    normalisation = normalisation_in;
  }

  std::string fprop_flav;
  std::string bprop_flav;
  std::string src_mom_prop;
  int gi;
  int gf;
  int gb;
  ::cvc::complex normalisation;

} oet_meson_twopt_meta_t;


/**
 * @brief Meta-description of a meson three-point function
 *
 * Describes the connected part of a three-point function
 *
 * \sum_{x_i, x_f, x_c}  e^{-i p_i x_i} e^{i p_f x_f} e^{i p_c x_c}
 * \bar{chi}_{f_1}(x_f) gf chi_{f_2}(x_f) \bar{chi}_{f_2}(x_c) gc chi_{f_2}(x_c)
 *     \bar{chi}_{f_2}(x_i) gi f_1(x_i)
 *
 * in the twisted basis. This is computed for fixed source-sink separation as
 *
 * Tr[ g5 (S_{f_1})^dag g5 gf g5 (S_{f_2})^dag g5 gc (S_{f_2}) gi ]
 * = -[ (D^{-1}_{f_2^dag} g5 gf g5 D^{-1}_{f_1^dag} gb eta ]^dag 
 *   g5 gc [D^{-1}_{f_2} gi eta] (1)
 *
 * using the property E(eta eta^dag) = 1 of the stochastic source field eta.
 *
 * In the present implementation, the two 'g5' bracketing 'gf' will be ignored 
 * in the final contraction, such that the behaviour of 'gf' under daggering and 
 * 'g5 gf g5' has to be included in the normalisation of the correlation function.
 *
 * As in the twopt function, gb may be just 'g5' or it might be more convenient
 * absorbe some gamma structure at the source location to minimize the number
 * of inversions required. Note that the g5 multiplying gf is implicitly included
 * in the contraction.
 *
 * In (1), the second factor is a propagator running from the source to the insertion
 * while the first factor is a sequential propagator: 
 * D^{-1}_{f_1^dag} gb eta runs from source, x_i, to sink, x_f, and 
 * (D^{-1}_{f_2^dag} g5 gf g5
 * from sink to insertion, x_c.
 *
 * Care has to be taken to use the correct phase factor at x_f because the
 * sequential propagator is daggered. Also, flavour change under gamma5 hermiticity
 * needs to be taken care of explicitly in the definition of the correlation function.
 *
 * As in the twopt function, the source momentum can be attached to the forward
 * or the backward propagator. In the latter case, the phase needs to be daggered.
 *
 * The base class carries the information on the backward and forward
 * propagator as well as the normalisation and the source and sink gamma
 * structures.
 *
 */
typedef struct oet_meson_threept_meta_t : oet_meson_twopt_meta_t 
{
  oet_meson_threept_meta_t(const std::string fprop_flav_in,
                     const std::string bprop_flav_in,
                     const std::string sprop_flav_in,
                     const std::string src_mom_prop_in,
                     const int gi_in,
                     const int gf_in,
                     const int gc_in,
                     const int gb_in,
                     const ::cvc::complex normalisation_in) :
    oet_meson_twopt_meta_t(fprop_flav_in, bprop_flav_in, src_mom_prop_in,
                           gi_in, gf_in, gb_in, normalisation_in )
  {
    gc = gc_in;
    sprop_flav = sprop_flav_in;
  }

  int gc;
  std::string sprop_flav;
} oet_meson_threept_meta_t; 

typedef struct MetaCollection {
  mom_lists_t mom_lists;

  std::shared_ptr< const std::vector<double> > ranspinor;

  int src_ts;
  std::map<std::string, stoch_prop_meta_t> props_meta;
  std::map<std::string, seq_stoch_prop_meta_t> seq_props_data;
  std::map<std::string, ts_stoch_src_meta_t> srcs_meta;
  DepGraph props_graph;
  
  DepGraph corrs_graph;
} MetaCollection;

typedef struct DataCollection {
  std::map<std::string, std::vector<double> > props_data;
  std::map<std::string, std::vector<double> > seq_props_data;
  std::map<std::string, std::vector<double> > deriv_props_data;
} DataCollection;

typedef struct OutputCollection {
  std::string corr_h5_filename;
  std::map< std::string, H5Correlator > corrs_data;
} OutputCollection;

//typedef struct threept_shifts_oet_meta_t : oet_meson_threept_meta_t
//{
//  threept_shifts_oet_meta_t(
//      const std::string fprop_flav_in,
//      const std::string bprop_flav_in,
//      const std::string sprop_flav_in,
//      const std::string src_mom_prop_in,
//      const int gi_in,
//      const int gf_in,
//      const int gc_in,
//      const int gb_in,
//      const std::vector<shift_t> left_shifts_in,
//      const std::vector<shift_t> right_shifts_in,
//      const ::cvc::complex normalisation_in
//      ) :
//    oet_meson_threept_meta_t(fprop_flav_in, bprop_flav_in, sprop_flav_in,
//                       src_mom_prop_in, gi_in, gf_in, gc_in, gb_in, normalisation_in)
//  {
//    left_shifts = left_shifts_in;
//    right_shifts = right_shifts_in; 
//  }
//  std::vector<shift_t> left_shifts;
//  std::vector<shift_t> right_shifts;
//} threept_shifts_oet_meta_t;

} // naemspace(cvc)
