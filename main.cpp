#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <random>


struct Process {
    std::string pid;
    int arrival_time;
    bool cpuBound;
    std::vector<int> cpu_bursts;
    std::vector<int> io_bursts;
};

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

std::vector<Process> generateProcesses(int n, int ncpu, double lambda, double upper_bound) {
    std::vector<Process> processes;

    for (int i = 0; i < n; i++) {
        Process p;
        p.cpuBound = i < ncpu;
        p.arrival_time = std::floor(next_exp(lambda, upper_bound));
        p.pid = std::string(1, 'A' + (i / 10)) + std::to_string(i % 10);

        int burstNum = std::ceil(drand48() * 32);
        std::cout << (p.cpuBound ? "CPU-bound" : "I/O-bound") << " process " << p.pid << ": arrival time " << p.arrival_time << "ms; " << burstNum << " CPU bursts:" << std::endl;

        for (int j = 0; j < burstNum; ++j) {
            int cpu_burst = std::ceil(next_exp(lambda, upper_bound));

            if (p.cpuBound) {
                cpu_burst *= 4;
            }

            p.cpu_bursts.push_back(cpu_burst);
            std::cout << "==> CPU burst " << cpu_burst << "ms";

            if (j < burstNum - 1) {  // don't make an io burst for last cpu burst
                int io_burst = std::ceil(next_exp(lambda, upper_bound)) * 8;
                
                if (p.cpuBound) {
                    io_burst /= 8;
                }

                p.io_bursts.push_back(io_burst);
                std::cout << " ==> I/O burst " << io_burst << "ms";
            }
            std::cout << std::endl;
        }
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

    std::vector<Process> Processes = generateProcesses(n, ncpu, lambda, upper_bound);

    return 0;
}