#include <iostream>
#include <vector>
#include <string>

class Process {
public:
    std::string pid;
    int arrival_time;
    bool cpuBound;
    std::vector<int> cpu_bursts;
    std::vector<int> io_bursts;

    Process(const std::string& pid, int arrival, bool cpuBound)
        : pid(pid), arrival_time(arrival), cpuBound(cpuBound) {}
};