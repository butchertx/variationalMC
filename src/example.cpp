/**
 * Example of a driver program, which can be used to run the pre-made example problems.
*/

#include <complex>
#include "vmctype.h"
#include "vmc_io.h"
#include "MemTimeTester.h"
// #include "Lattice.h"
// #include "MeanFieldAnsatz.h"
// #include "RandomEngine.h"
// #include "Wavefunction.h"
// #include "ProjectedState.h"
// #include "SpinModel.h"
// #include "VariationalMonteCarlo.h"
// #include "model_and_calculation_helper.h"

MemTimeTester timer;

int main(int argc, char* argv[]) {
    timer.flag_start_time("Total Program Time");
    std::string example_path, outfile_name;
    if (argc != 2) {
        std::cout << "Mandatory command line argument: <example directory>\n";
        std::cout << "Exiting...";
        return 1;
    }

    example_path = argv[1];
    std::cout << example_path << "\n";

    lattice_options lat_options = read_json_lattice_from_dir(example_path);
    mean_field_options wf_options = read_json_wavefunction_from_dir(example_path);
    model_options mdl_options = read_json_model_from_dir(example_path);
    vmc_options mc_options = read_json_vmc_from_dir(example_path);

    // read_json_full_input(&lat_options, &wf_options, &mdl_options, &mc_options, example_path);
    // makePath("./results");

    timer.flag_end_time("Total Program Time");
    timer.print_timers();

    return 0;

}