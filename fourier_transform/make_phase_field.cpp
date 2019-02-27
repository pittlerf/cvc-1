#include "global.h"
#include "loop_tools.h"
#include "cvc_complex.h"
#include "make_phase_field.hpp"
#include "types.h"

#include <vector>

namespace cvc {

  void make_phase_field(std::vector<::cvc::complex> & phases, const mom_t & p)
  {
    size_t const VOL3 = LX*LY*LZ;
    double const TWO_MPI = 2.0 * M_PI;
    double const px = TWO_MPI * (double)p.x / (double)LX_global;
    double const py = TWO_MPI * (double)p.y / (double)LY_global;
    double const pz = TWO_MPI * (double)p.z / (double)LZ_global;

    double const phase_offset = (double)( g_proc_coords[1] * LX ) * px +
                                (double)( g_proc_coords[2] * LY ) * py +
                                (double)( g_proc_coords[3] * LZ ) * pz;

    if( phases.size() != VOL3 ){
      phases.resize(VOL3);
    }
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
      double phase;
      unsigned int x, y, z;

      FOR_IN_PARALLEL(ix, 0, VOL3){
        x = g_lexic2coords[ix][1];
        y = g_lexic2coords[ix][2];
        z = g_lexic2coords[ix][3];

        phase = phase_offset + x*px + y*py + z*pz;
        phases[ix].re = cos(phase);
        phases[ix].im = sin(phase);
      }
    } // omp parallel region closing brace
  }

} // namespace(cvc)
