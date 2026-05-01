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

struct hnode;
typedef struct hnode hnode;

struct snode_mo;
typedef struct snode_mo snode_mo;

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
    
    // STMOA* specific fields
    unsigned tree_parent;
    hnode *hroot;              // Persistent leftist heap root for sidetracks
    bool hroot_computed;       // Lazy evaluation flag
    unsigned hroot_version;    // Version control for lazy heap cache invalidation
    
    // Perfect Heuristics / Pre-computation bounds
    unsigned g2_min_from_start; // Crystal Ball: absolute min g2 from start
    unsigned exact_h1;          // Perfect G1 heuristic from goal
    unsigned exact_h2;          // Perfect G2 heuristic from goal
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

// ---------------------------------------------------------
// Persistent Leftist Heap Node (represents a Sidetrack)
// ---------------------------------------------------------
struct hnode
{
    unsigned target_node;
    unsigned origin_node_id;  // ID of the node from which this sidetrack edge departs
    unsigned g1_delta;     // Detour penalty in G1
    int g2_delta;          // Detour penalty in G2
    
    hnode *left;
    hnode *right;
    int npl;               // Null Path Length (for leftist heap property)
    
    long long min_g2_term;      // Minimum possible G2 term for this sidetrack
    long long subtree_min_term; // Minimum possible G2 term in this entire subtree
};

// ---------------------------------------------------------
// Multi-Objective Search Node (used in Pareto Queue)
// ---------------------------------------------------------
struct snode_mo
{
    unsigned g1;               
    unsigned g2;                    
    double key;                
    unsigned long heapindex;   
    hnode *current_sidetrack;  // Pointer to the sidetrack applied to reach this node
};

#endif