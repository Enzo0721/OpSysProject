#ifndef PRCS_H
#define PRCS_H

#include <iostream>
#include <vector>
#include <string>

struct Process {
    std::string pid;
    int arrival_time;
    bool cpuBound;
    std::vector<int> cpu_bursts;
    std::vector<int> io_bursts;
};

#endif