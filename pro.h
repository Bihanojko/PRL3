// Project: Prirazeni poradi preorder vrcholum
// Author:  Nikola Valesova, xvales02
// Date:    25. 4. 2018

#include <mpi.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include <chrono>


typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;


#define TAG 0

// struct for storing MPI transaction result
MPI_Status stat;


// structure for storing information about one edge of given graph
struct EdgeRecord {
    int ID;             // sequence number of edge
    int vertices[2];    // ids of nodes which edge connects 
    int rever;          // index of edge reversed to current one
    int next;           // index of edge that is the next one starting from the same vertex
    bool isForward;     // true if edge is forward, false if it is retreat
};
