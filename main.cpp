#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <random>
#include <numeric>
#include <list>
#include <algorithm>

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

bool compare_by_arrival_time(const Process & left, const Process & right){
    return left.arrival_time < right.arrival_time;
}

void remove_latest(std::vector<Process> & start, std::vector<Process> & end){
    std::cout << "hello world" << std::endl;
}

void move_process(int time_change, std::vector<Process> & start, std::vector<Process> & destination){
                  //  ^^^^^^^^^^^^ for add tick + 1/2*tcs for context switches
    if (!start.empty()) {
        //int moving_time = tcs / 2;
        Process temp = start.front();
        start.erase(start.begin());

        temp.arrival_time = time_change;

        auto it = std::lower_bound(destination.begin(), destination.end(), temp, compare_by_arrival_time);
        destination.insert(it, temp);
    } else {
        std::cerr << "Error: Attempted to move a process from an empty vector.\n";
    }
}

std::string queue_string(const std::vector<Process> & queue){
    std::string full_string = "[Q";
    if(queue.size() > 0){
        for(int i = 0; i < queue.size(); ++i){
            full_string += " ";
            full_string += queue[i].pid;
        }
    }
    else{
        full_string += " empty";
    }

    full_string += "]";
    
    return full_string;
}

void remove_first_cpu_burst(std::vector<Process>& running_CPU_burst) {
    // Ensure that running_CPU_burst is not empty before trying to access its front
    if (!running_CPU_burst.empty()) {
        // Ensure that the process at the front has at least one CPU burst to remove
        if (!running_CPU_burst.front().cpu_bursts.empty()) {
            // Erase the first CPU burst
            running_CPU_burst.front().cpu_bursts.erase(running_CPU_burst.front().cpu_bursts.begin());
        } else {
            std::cerr << "Error: No CPU bursts to remove for the running process." << std::endl;
        }
    } else {
        std::cerr << "Error: No process is currently running." << std::endl;
    }
}

void remove_first_io_burst(std::vector<Process>& blocking_on_io) {
    // Ensure that blocking_on_io is not empty before trying to access its front
    if (!blocking_on_io.empty()) {
        // Ensure that the process at the front has at least one I/O burst to remove
        if (!blocking_on_io.front().io_bursts.empty()) {
            // Erase the first I/O burst
            blocking_on_io.front().io_bursts.erase(blocking_on_io.front().io_bursts.begin());
        } else {
            std::cerr << "Error: No I/O bursts to remove for the process." << std::endl;
        }
    } else {
        std::cerr << "Error: No process is currently blocking on I/O." << std::endl;
    }
}

void run_FCFS(std::vector<Process> & Processes, int tcs){
    int tick = 0;
    std::vector<Process> unarrived_processes;
    std::copy(Processes.begin(), Processes.end(), std::back_inserter(unarrived_processes));
    std::sort(unarrived_processes.begin(), unarrived_processes.end(), compare_by_arrival_time);

    int num_processes = Processes.size();

    std::vector<Process> queue;
    std::vector<Process> blocking_on_io;
    std::vector<Process> running_CPU_burst;

    std::vector<Process> completed; // To track completed processes

    int remaining_CS_t = 0;
    const int HALF_TCS = (int)tcs/2;

    std::cout << "time 0ms: Simulator started for FCFS [Q empty]" << std::endl;
    while(!unarrived_processes.empty() || !queue.empty() || !blocking_on_io.empty() || !running_CPU_burst.empty()){
        //To-do prevent processes from doing things while a context switch is occuring

        /*order:
        -Process arrival
        -Process starts using the CPU
        -Process finishes using the CPU
        -Process starts an I/O burst
        -Process finishes an I/O burst
        */
       std::string working_pid;

       //checking for newly arrived processes
       if(!unarrived_processes.empty()){
            while(!unarrived_processes.empty() && unarrived_processes.front().arrival_time == tick){
                working_pid = unarrived_processes.front().pid;
                move_process(tick + HALF_TCS, unarrived_processes, queue);                //semantically, the updated arrival time means the time it arrives in queue
                std::cout << "time " << tick << "ms: " << "Process "<< working_pid <<" arrived; added to ready queue " << queue_string(queue) <<std::endl;
            }
        }

        if (remaining_CS_t > 0) {
            remaining_CS_t--;
        }

        //moving from queue to CPU
        if(remaining_CS_t == 0 && queue.size() > 0 && running_CPU_burst.empty()){
            if(queue.front().cpu_bursts.size() != 0 && queue.front().arrival_time <= tick){
                working_pid = queue.front().pid; 
                int time_spent = queue.front().cpu_bursts.front();
                remaining_CS_t = HALF_TCS;

                //sets arrival time member to time that it completes burst
                move_process(tick + queue.front().cpu_bursts.front(), queue, running_CPU_burst);
                std::cout << "time " << tick << "ms: " << "Process "<< working_pid << " started using the CPU for " 
                << time_spent << "ms burst "  << queue_string(queue) << std::endl;

                remove_first_cpu_burst(running_CPU_burst);
            }
        }
        
        //check for if the CPU process is done
        if(!running_CPU_burst.empty() && running_CPU_burst[0].arrival_time == tick){
            std::cout <<"time " << tick << "ms: Process " <<running_CPU_burst[0].pid << " completed a CPU burst; " <<running_CPU_burst[0].cpu_bursts.size() <<" bursts to go " << queue_string(queue) <<std::endl;
        }
        //Move process out of CPU and into i/o or to completed list
        if(remaining_CS_t == 0 && !running_CPU_burst.empty() && running_CPU_burst[0].arrival_time <= tick){
            if(running_CPU_burst[0].cpu_bursts.empty()){
                // Process has completed all CPU bursts
                std::cout << "time " << tick << "ms: Process " << running_CPU_burst[0].pid << " terminated " << queue_string(queue) << std::endl;
                completed.push_back(running_CPU_burst.front());  // Move to completed list
                running_CPU_burst.clear(); // Clear the running CPU burst
            }
            else{
                std::cout <<"time " << tick << "ms: Process " << running_CPU_burst[0].pid
                <<" switching out of CPU; blocking on I/O until time " << tick + running_CPU_burst[0].io_bursts[0] + HALF_TCS << "ms " << queue_string(queue) <<std::endl;
                remaining_CS_t = tcs;
                move_process(tick + running_CPU_burst[0].io_bursts[0] + HALF_TCS, running_CPU_burst, blocking_on_io);
                //remove_first_io_burst(blocking_on_io); // Remove the first I/O burst
            }
                
        }
        
        //moving process out of io and back into queue
        if(!blocking_on_io.empty() && remaining_CS_t == 0 && blocking_on_io[0].arrival_time <= tick){
            remove_first_io_burst(blocking_on_io);
            working_pid = blocking_on_io[0].pid;
            move_process(tick + HALF_TCS, blocking_on_io, queue);
            std::cout <<"time " << tick << "ms: Process " << working_pid
            <<" completed I/O" << "; added to ready queue " << queue_string(queue) <<std::endl;
            remaining_CS_t = HALF_TCS;
        }

        tick++;
    }

    std::cout << "time " << tick << "ms: Simulator ended for FCFS" << std::endl;
}



int main(int argc, char** argv) {
    if (argc != 9) {
        std::cerr << "ERROR: Enter 7 command line arguments\n" << std::endl;
        return 1;
    }
    int n = std::stoi(argv[1]); /* Number of processes to simulate */
    int ncpu = std::stoi(argv[2]); /* Number of processes that are CPU-bound */
    int seed = std::stoi(argv[3]); /* seed for pseudo-random number sequence */
    double lambda = std::stod(argv[4]); /* lambda*/
    double upper_bound = std::stod(argv[5]); /* random value upper bound */
    int tcs = std::stoi(argv[6]); /* context switch time */
    float alpha = std::stof(argv[7]); /* true arrival time constant */
    int tslice = std::stoi(argv[8]); /* round robin time slice */

    if (n <= 0 || n > 260 || ncpu < 0 || ncpu > n || lambda <= 0 || upper_bound <= 0 
    || tcs <= 0 || (tcs % 2) != 0 || alpha > 1.0 || alpha < 0.0 || tslice <= 0) {
        std::cerr << "ERROR: Invalid input parameters\n" << std::endl;
        return 1;
    }
    srand48(seed);
    std::cout << "<<< PROJECT PART I" << std::endl;
    std::cout << "<<< -- process set (n=" << n << ") with " << ncpu << " CPU-bound process" << (ncpu == 1 ? "" : "es") << std::endl;
    std::cout << "<<< -- seed=" << seed << "; lambda=" << std::fixed << std::setprecision(6) << lambda << "; bound=" << std::setprecision(0) << upper_bound << std::endl;
    std::vector<Process> Processes = generateProcesses(n, ncpu, lambda, upper_bound);
    generateStatistics(Processes, n, ncpu);

    std::cout << "<<< PROJECT PART II"<< std::endl;
    std::cout << "<<< -- t_cs="<< tcs << "ms; alpha="<< std::fixed << std::setprecision(2) << alpha << "; t_slice=" << tslice <<"ms" << std::endl;
    run_FCFS(Processes, tcs);
    return 0;
}