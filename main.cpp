#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "ERROR: Enter 5 command line arguments\n";
        return 1;
    }

    int n = std::stoi(argv[1]);
    int ncpu = std::stoi(argv[2]);
    int seed = std::stoi(argv[3]);
    double lambda = std::stod(argv[4]);
    double upper_bound = std::stod(argv[5]);

    if (n <= 0 || ncpu < 0 || ncpu > n || lambda <= 0 || upper_bound <= 0) {
        std::cerr << "ERROR: Invalid input parameters\n";
        return 1;
    }

    std::cout << "<<< PROJECT PART I" << std::endl;
    std::cout << "<<< -- process set(n=" << n << ") with " << ncpu << " CPU-bound process" << (ncpu == 1 ? "" : "es") << std::endl;
    std::cout << "<<< -- seed=" << seed << "; lambda=" << std::fixed << std::setprecision(6) << lambda << "; bound=" << std::setprecision(0) << upper_bound << std::endl;

}
