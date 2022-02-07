#include <complex>
#include "../src/vmctype.h"
#include "../src/vmc_io.h"
#include "../src/MemTimeTester.h"
#include "../src/Lattice.h"
#include "../src/RandomEngine.h"
#include "../src/Wavefunction.h"

using namespace vmctype;
MemTimeTester timer;

bool test_jastrow_sz();
bool test_jastrow_sz2();

int main(int argc, char* argv[]) {
    timer.flag_start_time("Total Test Time");

    test_jastrow_sz();

    timer.flag_end_time("Total Test Time");
}

bool test_jastrow_sz() {
    bool passed = true;

    int num_tests = 1000;
    JastrowFactorOptions j;
    j.distance_max = 1;
    j.values = { -0.1, 0.2 };

    Lattice lattice(Lattice_type_from_string("triangle"), vec3<int>(6, 6, 1), vec3<int>(1, 1, 1));
    RandomEngine rand(-1, lattice.get_N(), 1);
    int N = lattice.get_N(), N0 = lattice.get_N() / 3;
    std::vector<int> configuration = rand.get_rand_spin_state(std::vector<int>{ (N - N0) / 2, N0, (N - N0) / 2 }, N);

    std::vector<JastrowFactor> jlist;
    std::vector<std::vector<int>> neighbors;
    for (int i = 0; i < j.distance_max; ++i) {
        neighbors = lattice.get_neighbors(i);
        jlist.push_back(JastrowFactor(j.values[i], neighbors));
    }

    JastrowTable jtable(jlist);
    jtable.initialize_tables(configuration);

    std::vector<int> sites2(2), spins2(2), sites3(3), spins3(3);
    double oldval, newval, lazyratio;

    //test lazy eval
    for (int i = 0; i < num_tests; ++i) {
        oldval = jtable.greedy_eval(configuration);
        sites2 = { rand.get_rand_site(), rand.get_rand_site() };
        spins2 = { configuration[sites2[1]], configuration[sites2[0]] };
        lazyratio = jtable.lazy_eval(sites2, spins2, configuration);
        configuration[sites2[0]] = spins2[0];
        configuration[sites2[1]] = spins2[1];
        newval = jtable.greedy_eval(configuration);
        if (std::abs(lazyratio - newval / oldval) > 1e-6) {
            std::cout << "spins = " << vec2str(spins2) << "\n";
            std::cout << "Lazy = " << lazyratio << ", greedy = " << newval / oldval << "\n";
            passed = false;
        }
    }

    //test update
    jtable.initialize_tables(configuration);
    for (int i = 0; i < num_tests; ++i) {
        oldval = jtable.greedy_eval(configuration);
        sites2 = { rand.get_rand_site(), rand.get_rand_site() };
        spins2 = { configuration[sites2[1]], configuration[sites2[0]] };
        lazyratio = jtable.lazy_eval(sites2, spins2, configuration);
        configuration[sites2[0]] = spins2[0];
        configuration[sites2[1]] = spins2[1];
        newval = jtable.greedy_eval(configuration);
        if (std::abs(lazyratio - newval / oldval) > 1e-6) {
            std::cout << "spins = " << vec2str(spins2) << "\n";
            std::cout << "Lazy = " << lazyratio << ", greedy = " << newval / oldval << "\n";
        }
    }

    return passed;
}