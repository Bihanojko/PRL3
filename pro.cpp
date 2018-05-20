// Project: Prirazeni poradi preorder vrcholum
// Author:  Nikola Valesova, xvales02
// Date:    25. 4. 2018


#include "pro.h"

using namespace std;


// create new edge and initialize it with given values
EdgeRecord initializeEdge(int ID, int from, int to, int rever, int next, bool isForward) {
    EdgeRecord Edge;
    Edge.ID = ID;
    Edge.vertices[0] = from;
    Edge.vertices[1] = to;
    Edge.rever = rever;
    Edge.next = next;
    Edge.isForward = isForward;
    return Edge;
}

// create and return an initialized new edge 
EdgeRecord createEdge(int type, int nodeCount, int graphLength, int i) {
    // edge to parent node
    if (type == 0) {
        if (2 * i + 1 < nodeCount)
            return initializeEdge(i + nodeCount - 2, i, (i - 1) / 2, i - 1, graphLength + 1, false);
        else
            return initializeEdge(i + nodeCount - 2, i, (i - 1) / 2, i - 1, -1, false);
    }
    // edge to left son
    else if (type == 1) {
        if (i * 2 + 2 < nodeCount)
            return initializeEdge(i * 2, i, i * 2 + 1, i * 2 + nodeCount - 1, graphLength + 1, true);
        else
            return initializeEdge(i * 2, i, i * 2 + 1, i * 2 + nodeCount - 1, -1, true);
    }
    // edge to right son
    else if (type == 2)
        return initializeEdge(i * 2 + 1, i, i * 2 + 2, i * 2 + nodeCount, -1, true);
}

// read and process input sequence and create adjacency list and data structure for storing the graph
void createGraph(int nodeCount, vector<EdgeRecord> &graph, int *AdjList) {
    for (int i = 0; i < nodeCount; i++) {
        AdjList[i] = graph.size();

        // if current node has a parent, append edge to parent node
        if (i != 0)
            graph.push_back(createEdge(0, nodeCount, graph.size(), i));
        // if current node has a left son, append edge to left son node
        if (i * 2 + 1 < nodeCount)
            graph.push_back(createEdge(1, nodeCount, graph.size(), i));
        // if current node has a right son, append edge to right son node
        if (i * 2 + 2 < nodeCount)
            graph.push_back(createEdge(2, nodeCount, graph.size(), i));
    }
}

// set correct IDs of reverse edges
void setReverseEdges(vector<EdgeRecord> &graph) {
    int *correctedRevs = new int [graph.size()];

    for (int i = 0; i < graph.size(); i++) {
        for (int j = 0; j < graph.size(); j++) {
            if (graph[i].rever == graph[j].ID) {
                correctedRevs[i] = j;
                break;
            }
        }
    }

    for (int i = 0; i < graph.size(); i++)
        graph[i].rever = correctedRevs[i];

    delete[] correctedRevs;
}

// receive one value from processor with ID neighID 
int receiveValue(int neighID) {
    int value;
    MPI_Recv(&value, 1, MPI_INT, neighID, TAG, MPI_COMM_WORLD, &stat);
    return value;
}

// send one value to processor with ID neighID
void sendValue(int neighID, int value) {
    MPI_Send(&value, 1, MPI_INT, neighID, TAG, MPI_COMM_WORLD);
}

// receive size of data being transferred and receive the data itself
void receiveData(int neighID, int *data) {
    int size;
    MPI_Recv(&size, 1, MPI_INT, neighID, TAG, MPI_COMM_WORLD, &stat);
    MPI_Recv(data, size, MPI_INT, neighID, TAG, MPI_COMM_WORLD, &stat);
}

// send size of data being transferred and send the data itself
void sendData(int neighID, int numprocs, int *data) {
    MPI_Send(&numprocs, 1, MPI_INT, neighID, TAG, MPI_COMM_WORLD);
    MPI_Send(data, numprocs, MPI_INT, neighID, TAG, MPI_COMM_WORLD);
}

// receive values from all processors and fill whole data array, then execute broadcast algorithm
void receiveAndBroadcast(int myID, int numprocs, int value, int *data) {
    data[myID] = value;

    // receive information to fill values into data array
    for (int i = 1; i < numprocs; i++)
        data[i] = receiveValue(i);

    // broadcast algorithm to simulate simultaneous reading
    for (int i = 1; i < numprocs; i *= 2)
        sendData(i, numprocs, data);
}

// send my value to processors with ID 0 and receive the filled data array, then execute broadcast algorithm
void sendAndBroadcast(int myID, int numprocs, int value, int *data) {
    // send my value to processor with ID 0
    sendValue(0, value);
    
    // receive filled data array using broadcast algorithm
    receiveData(myID - pow(2, trunc(log2(myID))), data);

    // continue in broadcasting the filled data array to next neighbour processors
    for (int i = pow(2, floor(log2(myID)) + 1); i + myID < numprocs; i *= 2)
        sendData(i + myID, numprocs, data);
}

// fill data array with values from all processors and then broadcast the whole array to all processors
void exchangeInformation(int myID, int numprocs, int value, int *data) {    
    if (myID == 0)
        receiveAndBroadcast(myID, numprocs, value, data);
    else
        sendAndBroadcast(myID, numprocs, value, data);
}

// computation of suffix sum and correction of edge order
void suffixSum(int myID, int numprocs, bool isForward, int *Rank, int *EtourLink) {
    int iterationCount = trunc(log2(numprocs)) + 1;
    int weight;
    int Etour;

    // set initial weight of edge
    if (isForward)
        weight = 1;
    else
        weight = 0;

    for (int i = 0; i < iterationCount; i++) {
        exchangeInformation(myID, numprocs, weight, Rank);
        weight += Rank[EtourLink[myID]];
        Etour = EtourLink[EtourLink[myID]];
        exchangeInformation(myID, numprocs, Etour, EtourLink);
    }

    weight = numprocs / 2 - weight + 1;
    exchangeInformation(myID, numprocs, weight, Rank);
}

// MAIN
int main(int argc, char *argv[]) {
    int numprocs;                       // number of processors
    int myID;                           // my processor id
    int nodeCount = strlen(argv[1]);    // number of nodes of given graph

    // MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myID);

    vector<EdgeRecord> graph;                       // vector of edges, inner representation of graph
    int *AdjList = new int [nodeCount];             // adjacency list for Euler tour
    int Etour;                                      // etour(e) of current edge
    int *EtourLink = new int [numprocs];            // array to store Euler Tour of all processors
    int *Rank = new int [numprocs];                 // array of ranks of Euler Tour algorithm
    int *sortedSequence = new int [nodeCount];    // array to store values of nodes during preorder traversal 

    // graph contains only one node, no one to exchange information with
    if (nodeCount == 1)
        cout << argv[1][0] << endl;
    // otherwise execute preorder traversal
    else {
        // create adjacency list and inner graph representation
        createGraph(nodeCount, graph, AdjList);
        setReverseEdges(graph);

        // if edge reversed to current edge has a 'next'   
        if (graph[graph[myID].rever].next != -1)
            // set Etour(e) to next(e_R)
            Etour = graph[graph[myID].rever].next;
        else
            // otherwise set Etour(e) to to the first item of AdjList(v)
            Etour = AdjList[graph[myID].vertices[1]];
        
        // if current edge is edge going to the root node, set its Etour(e) to e
        if (Etour == 0)
            Etour = myID;

        // construct whole Euler Tour array
        exchangeInformation(myID, numprocs, Etour, EtourLink);

        // compute suffix sum
        suffixSum(myID, numprocs, graph[myID].isForward, Rank, EtourLink);

        // create ordered sequence and output it
        if (myID == 0) {
            sortedSequence[Rank[myID]] = (int)argv[1][graph[myID].vertices[1]];
            sortedSequence[0] = (int)argv[1][0];
            sendData(1, nodeCount, sortedSequence);
            receiveData(numprocs - 1, sortedSequence);

            // output values in right order
            if (myID == 0) {
                for (int i = 0; i < nodeCount; i++)
                    cout << (char)sortedSequence[i];
                cout << endl;
            }
        }
        // receive data, add your value to the right position and pass on
        else {
            receiveData(myID - 1, sortedSequence);            
            if (graph[myID].isForward)
                sortedSequence[Rank[myID]] = (int)argv[1][graph[myID].vertices[1]];
            if (myID + 1 < numprocs)
                sendData(myID + 1, nodeCount, sortedSequence);
            else
                sendData(0, nodeCount, sortedSequence);
        }
    }

    delete[] AdjList;
    delete[] EtourLink;
    delete[] Rank;
    delete[] sortedSequence;
    MPI_Finalize();
    return 0;
}
