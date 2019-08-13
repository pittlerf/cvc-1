#ifndef ENUMS_HPP
#define ENUMS_HPP

namespace cvc {

typedef enum latDim_t {
  DIM_T = 0,
  DIM_X,
  DIM_Y,
  DIM_Z,
  DIM_NDIM
} latDim_t;

extern const char latDim_names[];

typedef enum shift_dir_t {
  DIR_FWD = 0,
  DIR_BWD,
  DIR_NDIR
} shift_dir_t;

extern const char shift_dir_names[];

typedef enum ExitCode_t {
  CVC_EXIT_SUCCESS = 0,
  CVC_EXIT_CORE_INIT_FAILURE = 1,
  CVC_EXIT_H5_KEY_LENGTH = 2,
  CVC_EXIT_GAUGE_INIT_FAILURE = 3,
  CVC_EXIT_UTIL_FUNCTION_FAILURE = 4,
  CVC_EXIT_MALLOC_FAILURE = 5,
  CVC_EXIT_MPI_FAILURE = 6,
  CVC_EXIT_INVALID_INPUT = 7,
  CVC_EXIT_SNPRINTF_OVERFLOW = 8,
  CVC_EXIT_NCODES
} ExitCode_t;

typedef enum NoiseType_t {
  GAUSSIAN_NOISE = 1,
  Z2_NOISE,
  N_NOISE_TYPES
} NoiseType_t;

typedef enum QuarkSmearingType_t {
  QUARK_SMEAR_NONE = 0,
  QUARK_SMEAR_JACOBI,
  QUARK_SMEAR_MOM_JACOBI,
  QUARK_SMEAR_NTYPES
} QuarkSmearingType_t;

typedef enum GaugeSmearingType_t {
  GAUGE_SMEAR_NONE = 0,
  GAUGE_SMEAR_APE,
  GAUGE_SMEAR_NTYPES
} GaugeSmearingType_t;

typedef enum cvc_GammaBasis {
  UKQCD = 0,
  CHIRAL_TMLQCD,
} cvc_GammaBasis;
  
}

#endif
