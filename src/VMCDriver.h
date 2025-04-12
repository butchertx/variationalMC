#pragma once
#include "vmctype.h"

using namespace vmctype;

class VMCResults;

class VMCDriver {

    LatticeOptions lat_options;
    WavefunctionOptions wf_options;
    ModelOptions mdl_options;
    VMCOptions mc_options;

public:

    VMCDriver(LatticeOptions lat_opts, WavefunctionOptions wf_opts, ModelOptions mdl_opts, VMCOptions mc_opts)
        : lat_options(lat_opts), wf_options(wf_opts), mdl_options(mdl_opts), mc_options(mc_opts) {}

    VMCResults run();

};