#ifndef _PDR_APPROXIMATOR_HPP
#define _PDR_APPROXIMATOR_HPP

#include "DummyApproximator.hpp"
#ifdef PRODM_HAVE_HPEZ
#include "HPEZApproximator.hpp"
#include "GEApproximator.hpp"
#endif
#ifdef PRODM_HAVE_SZ3
#include "SZ3Approximator.hpp"
#endif
#ifdef PRODM_HAVE_SZ2
#include "SZ2Approximator.hpp"
#endif
#ifdef PRODM_HAVE_MGARD
#include "MGARDApproximator.hpp"
#endif

#endif
