#ifndef SG_BOASTAR_H
#define SG_BOASTAR_H

#define MAX_SOLUTIONS 1000000
#define MAX_RECYCLE   1000000

extern gnode *graph_node;
extern unsigned num_gnodes;
extern unsigned adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned pred_adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned goal, start; 
extern unsigned solutions[MAX_SOLUTIONS][2];
extern unsigned nsolutions;
extern unsigned long long int stat_generated;
extern unsigned long long int stat_expansions;
extern unsigned long long int stat_exsarcted;
extern unsigned long long int stat_created;
extern unsigned long long int stat_recycled;
extern unsigned stat_pruned_sum;
extern unsigned stat_pruned;

void initialize_parameters();
int backward_dijkstra(int dim);
int sg_boastar();
void call_sg_boastar();

#endif
