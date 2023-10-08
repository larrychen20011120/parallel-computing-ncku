#include <iostream>
#include <fstream>
#include <string>
#include "mpi.h"

#define MAX_BIT_SIZE 32
#define MAIN_PROCESS 0
#define CODE 0

using namespace std;

int main(int argc, char *argv[]) {

    // problem variables
    int n, m;
    int count, cost, part;
    unsigned int allChoices = 1, endCase = 1;
    unsigned int localCount = 0;
    string filename;
    // allocate the memory of all test task
    unsigned int tasks[MAX_BIT_SIZE] = {0};

    // mpi variables
    int currentId;
    int processCount;

    // start mpi
    MPI_Init(&argc, &argv);
    // retrieve id
    MPI_Comm_rank(MPI_COMM_WORLD, &currentId);
    // get the total process count
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    // running command contains input file and output file
    // use cin/cout can handle the standard I/O from running command
    if (currentId == MAIN_PROCESS) {
        // readfile
        cin >> filename;
        ifstream ifs(filename);

        // read n, m values
        ifs >> n >> m;

        // 2^m possible choices
        for (int i = 0; i < m; i++)
            allChoices *= 2;
        for (int i = 0; i < n; i++)
            endCase *= 2;
        // convert to all bit 1
        --endCase;

        for (int i = 0; i < m; i++) {
            ifs >> count >> cost;
            for (int j = 0; j < count; j++) {
                ifs >> part;
                // set part start from 0
                --part;
                // the part-th set to 1
                tasks[i] |= ((unsigned int)1 << part);
            }
        }
        ifs.close();
    }

    // broadcasting necessary data
    MPI_Bcast(&m, 1, MPI_INT, MAIN_PROCESS, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, MAIN_PROCESS, MPI_COMM_WORLD);
    MPI_Bcast(&allChoices, 1, MPI_INT, MAIN_PROCESS, MPI_COMM_WORLD);
    MPI_Bcast(&endCase, 1, MPI_INT, MAIN_PROCESS, MPI_COMM_WORLD);
    MPI_Bcast(&tasks, MAX_BIT_SIZE, MPI_UNSIGNED, 
              MAIN_PROCESS, MPI_COMM_WORLD);
    
    // synchronize
    MPI_Barrier(MPI_COMM_WORLD);

    // divide all choices to all processes
    for (unsigned int currentChoice = currentId; currentChoice < allChoices;
         currentChoice += processCount) {

        // myset will join all sets together where current bit is 1
        unsigned int myset = 0;

        for (int j = 0; j < m; j++) {
            if (currentChoice & ((unsigned int)1 << j)) {
                // if the j-th bit is 1,
                // then join myset with that set
                myset |= tasks[j];
            }
        }
        if (myset == endCase) // this means contains all n partitions
            localCount++;
    }
    
    // synchronize
    MPI_Barrier(MPI_COMM_WORLD);
    
    // output result
    if (currentId != MAIN_PROCESS) { // sender
        MPI_Send(&localCount, 1, MPI_UNSIGNED, MAIN_PROCESS,
                 CODE, MPI_COMM_WORLD);
    } else { // receiver
        for (int id = 1; id < processCount; id++) {
            unsigned int recvData;
            MPI_Recv(&recvData, 1, MPI_UNSIGNED, id, 
                     CODE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            localCount += recvData;
        } 
    }

    if (currentId == MAIN_PROCESS)
        cout << localCount;

    // end mpi
    MPI_Finalize();
    
    return 0;
}