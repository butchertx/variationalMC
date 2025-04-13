#include <gtest/gtest.h>
#include <filesystem>
#include "mkl_types.h"

#include <vmctype.h>
#include <vmc_io.h>
#include <Lattice.h>
#include <VMCDriver.h>
#include <VMCResults.h>

// Global Variables

#ifdef EXAMPLES_PATH
    const std::string EXAMPLES_DIR = EXAMPLES_PATH;
#else
    const std::string EXAMPLES_DIR = std::filesystem::absolute("../../examples/");
#endif

TEST(FullExampleEnergyTest, SpinHalf1dTrivialExample) {
    std::string example_path = EXAMPLES_DIR + "spin_half/1d/trivial";
    vmctype::LatticeOptions lat_options = read_json_lattice_from_dir(example_path);
    WavefunctionOptions wf_options = read_json_wavefunction_from_dir(example_path);
    ModelOptions mdl_options = read_json_model_from_dir(example_path);
    VMCOptions mc_options = read_json_vmc_from_dir(example_path);
    VMCDriver driver(lat_options, wf_options, mdl_options, mc_options);
    VMCResults results = driver.run();
    std::complex<double> E = results.get_energy();
    std::complex<double> E_err = results.get_energy_err();

    std::cout << "Energy = " << E << " +/- " << E_err << "\n";
    EXPECT_NEAR(std::real(E), 0.375, 0.01);
}

TEST(FullExampleEnergyTest, SpinOne1dTrivialExample) {
    std::string example_path = EXAMPLES_DIR + "spin_one/1d/trivial";
    vmctype::LatticeOptions lat_options = read_json_lattice_from_dir(example_path);
    WavefunctionOptions wf_options = read_json_wavefunction_from_dir(example_path);
    ModelOptions mdl_options = read_json_model_from_dir(example_path);
    VMCOptions mc_options = read_json_vmc_from_dir(example_path);
    VMCDriver driver(lat_options, wf_options, mdl_options, mc_options);
    VMCResults results = driver.run();
    std::complex<double> E = results.get_energy();
    std::complex<double> E_err = results.get_energy_err();

    std::cout << "Energy = " << E << " +/- " << E_err << "\n";
    EXPECT_NEAR(std::real(E), -1.0, 0.01);
}