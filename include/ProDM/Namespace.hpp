#include "ProDM/Utils/StatUtils.hpp"
#ifndef PRODM_NAMESPACE_HPP
#define PRODM_NAMESPACE_HPP

// ProDM namespace layout.
//
// Every component of the library lives under the umbrella namespace ProDM:
//   ProDM         the machinery shared by both pipelines: bitplane encoders, level compressors, error control
//                 (interfaces, collectors, the linear estimator, size interpreters), retrievers, writers, utilities
//   ProDM::MDR    multilevel (SC'21 / SC'26) pipeline: decomposers, interleavers, tuners, refactors, reconstructors,
//                 and the estimators whose constants come from the multilevel bases
//   ProDM::PDR    approximation-based (TVCG'23 / QProR) pipeline: approximators, refactors, reconstructors
//   ProDM::MGARDx the in-house multilevel decomposition internals (Decomposer/MultiLevel/MGARDx)
//   ProDM::readfile, writefile, print_statistics  library-wide file and statistics helpers (Utils/IOUtils.hpp)
//   ProDM::Legacy code kept only to reproduce prior papers (GE synthesizer recipes, WeightReconstructor, QoIRefactor)
//
// The unqualified names MDR and PDR used by existing applications remain valid through the
// aliases below. Headers reopen the namespaces with the nested form (namespace ProDM::MDR { ... });
// a header must include this file before that, so that the aliases are visible everywhere.
namespace ProDM {
    namespace MDR {}
    namespace PDR {}
    namespace MGARDx {}
    namespace Legacy {}
}
namespace MDR = ProDM::MDR;
namespace PDR = ProDM::PDR;

#endif
