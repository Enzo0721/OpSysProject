#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <random>
#include <numeric>

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
        std::cout << (p.cpuBound ? "CPU-bound" : "I/O-bound") << " process " << p.pid 
                << ": arrival time " << p.arrival_time << "ms; " << burstNum << " CPU " << (burstNum == 1 ? "burst" : "bursts") << ":" << std::endl;

        for (int j = 0; j < burstNum; j++) {
            int cpu_burst = std::ceil(next_exp(lambda, upper_bound));
            if (p.cpuBound) cpu_burst *= 4;
            p.cpu_bursts.push_back(cpu_burst);
            std::cout << "==> CPU burst " << cpu_burst << "ms";
            if (j < burstNum - 1) {  // don't make an io burst for last cpu burst
                int io_burst = std::ceil(next_exp(lambda, upper_bound)) * 8;
                if (p.cpuBound) io_burst /= 8;
                p.io_bursts.push_back(io_burst);
                std::cout << " ==> I/O burst " << io_burst << "ms";
            }
            std::cout << std::endl;
        }
        processes.push_back(p);
    }
    return processes;
}

void set_nan_zero(double & val){
    val = std::isnan(val) ? 0.0 : val;
}

void generateStatistics(const std::vector<Process>& processes, int n, int nCpu) {
    std::ofstream outFile("simout.txt");
    outFile << std::fixed << std::setprecision(3);
    int nIo = n - nCpu;
    double cpuBoundCpuTotal = 0.0, ioBoundCpuTotal = 0.0, cpuBoundIoTotal = 0.0, ioBoundIoTotal = 0.0;
    int cpuBoundCpuCount = 0, ioBoundCpuCount = 0, cpuBoundIoCount = 0, ioBoundIoCount = 0;

    for (const auto& p : processes) {
        if (p.cpuBound) {
            cpuBoundCpuTotal += std::accumulate(p.cpu_bursts.begin(), p.cpu_bursts.end(), 0.0);
            cpuBoundIoTotal += std::accumulate(p.io_bursts.begin(), p.io_bursts.end(), 0.0);
            cpuBoundCpuCount += p.cpu_bursts.size();
            cpuBoundIoCount += p.io_bursts.size();
        } else {
            ioBoundCpuTotal += std::accumulate(p.cpu_bursts.begin(), p.cpu_bursts.end(), 0.0);
            ioBoundIoTotal += std::accumulate(p.io_bursts.begin(), p.io_bursts.end(), 0.0);
            ioBoundCpuCount += p.cpu_bursts.size();
            ioBoundIoCount += p.io_bursts.size();
        }
    }
    
    double cpuBoundAvgCpu = std::ceil(cpuBoundCpuTotal / cpuBoundCpuCount * 1000) / 1000;
    set_nan_zero(cpuBoundAvgCpu);
    double ioBoundAvgCpu = std::ceil(ioBoundCpuTotal / ioBoundCpuCount * 1000) / 1000;
    set_nan_zero(ioBoundAvgCpu);
    double overallAvgCpu = std::ceil((cpuBoundCpuTotal + ioBoundCpuTotal) / (cpuBoundCpuCount + ioBoundCpuCount) * 1000) / 1000;

    double cpuBoundAvgIo = std::ceil(cpuBoundIoTotal / cpuBoundIoCount * 1000) / 1000;
    set_nan_zero(cpuBoundAvgIo);
    double ioBoundAvgIo = std::ceil(ioBoundIoTotal / ioBoundIoCount * 1000) / 1000;
    set_nan_zero(ioBoundAvgIo);
    double overallAvgIo = std::ceil((cpuBoundIoTotal + ioBoundIoTotal) / (cpuBoundIoCount + ioBoundIoCount) * 1000) / 1000;

    outFile << "-- number of processes: " << n << std::endl;
    outFile << "-- number of CPU-bound processes: " << nCpu << std::endl;
    outFile << "-- number of I/O-bound processes: " << nIo << std::endl;
    outFile << "-- CPU-bound average CPU burst time: " << cpuBoundAvgCpu << " ms" << std::endl;
    outFile << "-- I/O-bound average CPU burst time: " << ioBoundAvgCpu << " ms" << std::endl;
    outFile << "-- overall average CPU burst time: " << overallAvgCpu << " ms" << std::endl;
    outFile << "-- CPU-bound average I/O burst time: " << cpuBoundAvgIo << " ms" << std::endl;
    outFile << "-- I/O-bound average I/O burst time: " << ioBoundAvgIo << " ms" << std::endl;
    outFile << "-- overall average I/O burst time: " << overallAvgIo << " ms" << std::endl;

    outFile.close();
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "ERROR: Enter 5 command line arguments\n" << std::endl;
        return 1;
    }
    int n = std::stoi(argv[1]); /* Number of processes to simulate */
    int ncpu = std::stoi(argv[2]); /* Number of processes that are CPU-bound */
    int seed = std::stoi(argv[3]); /* seed for pseudo-random number sequence */
    double lambda = std::stod(argv[4]); /* lambda*/
    double upper_bound = std::stod(argv[5]); /* random value upper bound */
    if (n <= 0 || n > 260 || ncpu < 0 || ncpu > n || lambda <= 0 || upper_bound <= 0) {
        std::cerr << "ERROR: Invalid input parameters\n" << std::endl;
        return 1;
    }
    srand48(seed);
    std::cout << "<<< PROJECT PART I" << std::endl;
    std::cout << "<<< -- process set (n=" << n << ") with " << ncpu << " CPU-bound process" << (ncpu == 1 ? "" : "es") << std::endl;
    std::cout << "<<< -- seed=" << seed << "; lambda=" << std::fixed << std::setprecision(6) << lambda << "; bound=" << std::setprecision(0) << upper_bound << std::endl;
    std::vector<Process> Processes = generateProcesses(n, ncpu, lambda, upper_bound);
    generateStatistics(Processes, n, ncpu);
    return 0;
}