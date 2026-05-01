#ifndef SG_BOASTAR_H
#define SG_BOASTAR_H

#include "node.h"

extern gnode *graph_node;
extern unsigned num_gnodes;
extern unsigned adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned pred_adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned goal, start;


#define MAX_SOLUTIONS 1000000
#define MAX_RECYCLE   1000000
#define MAX_WARM_STARTS 5

extern unsigned solutions[MAX_SOLUTIONS][2];
extern unsigned nsolutions;
extern unsigned long long int stat_generated;
extern unsigned long long int stat_expansions;
extern unsigned long long int stat_exsarcted;
extern unsigned long long int stat_created;
extern unsigned long long int stat_recycled;
extern unsigned stat_pruned;

// Warm-start tracking
extern unsigned ws_g1[MAX_WARM_STARTS];
extern unsigned ws_g2[MAX_WARM_STARTS];
extern int warm_start_count;

void initialize_parameters();
int backward_dijkstra(int dim);
void find_warm_start_solutions();
int sg_boastar();
void call_sg_boastar();

#endif