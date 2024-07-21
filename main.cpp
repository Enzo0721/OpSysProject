#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <random>
#include "process.h"

double next_exp(double lambda, double upper_bound) {
    double exp_val;
    do {
        double u;
        do {
            u = drand48();
        } while (u == 0); // To avoid log(0)
        exp_val = -log(u) / lambda;
    } while (exp_val > upper_bound);
    return exp_val;
}

std::vector<Process>  generateProcess(int n, int ncpu, double lambda, double upper_bound) {
    std::vector <Process> processes;

    for (int i = 0; i < n; i++) {
        Process p;
        p.pid = std::string(1, 'A' + i / 10) + std::to_string(i % 10);
        p.cpuBound = i < ncpu;
        p.arrival_time = std::floor(next_exp(lambda, upper_bound));
        std::cout << "PROCESS ID: " << p.pid << std::endl;
        std::cout << "PROCESS ARRIVAL: " << p.arrival_time << std::endl;
        int burstNum = std::ceil(drand48() * 32);
        std::cout << "PROCESS BURST AMOUNT: " << burstNum << std::endl;
        processes.push_back(p);
    } 
    return processes;
}   

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "ERROR: Enter 5 command line arguments\n";
        return 1;
    }

    int n = std::stoi(argv[1]); /* Number of processes to simulate */
    int ncpu = std::stoi(argv[2]); /* Number of processes that are CPU-bound */
    int seed = std::stoi(argv[3]); /* seed for pseudo-random number sequence */
    double lambda = std::stod(argv[4]); /* lambda*/
    double upper_bound = std::stod(argv[5]); /* random value upper bound */

    if (n <= 0 || ncpu < 0 || ncpu > n || lambda <= 0 || upper_bound <= 0) {
        std::cerr << "ERROR: Invalid input parameters\n";
        return 1;
    }

    srand48(seed);

    std::cout << "<<< PROJECT PART I" << std::endl;
    std::cout << "<<< -- process set(n=" << n << ") with " << ncpu << " CPU-bound process" << (ncpu == 1 ? "" : "es") << std::endl;
    std::cout << "<<< -- seed=" << seed << "; lambda=" << std::fixed << std::setprecision(6) << lambda << "; bound=" << std::setprecision(0) << upper_bound << std::endl;

    std::vector<Process> Processes = generateProcess(n, ncpu, lambda, upper_bound);

    return 0;
}
