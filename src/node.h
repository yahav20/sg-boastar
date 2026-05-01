/////////////////////////////////////////////////////////////////////
// node.h
// Structure Definitions for Search Algorithms
/////////////////////////////////////////////////////////////////////

#ifndef MAZEH
#define MAZEH

#include "include.h"
#include <stdbool.h>

struct gnode;
typedef struct gnode gnode;

struct snode;
typedef struct snode snode;

// ---------------------------------------------------------
// Graph Node Representation
// ---------------------------------------------------------
struct gnode 
{
    long long int id;
    
    // Exact objective costs from start
    unsigned g1;
    unsigned g2;
    
    // Base heuristic info (if needed by other algorithms)
    unsigned h1;
    unsigned h2;
    double lat;
    double lon;
    
    // A* and Dijkstra internals
    unsigned long long int key;
    unsigned gmin;             // Used for suffix domination pruning
    unsigned long heapindex;
 
    snode *gopfirst;
    snode *goplast;
};

// ---------------------------------------------------------
// Standard Search Node (used by BOA*)
// ---------------------------------------------------------
struct snode 
{
    int state;
    unsigned g1;
    unsigned g2;
    double key;
    unsigned long heapindex;
    snode *searchtree;
    snode *gopnext;
    snode *gopprev;
};


#endif