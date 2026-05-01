#include "heap.h"
#include "node.h"
#include "include.h"
#include "sg_boastar.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdbool.h>
#define MAX_SNODES 10000000

extern bool is_timeout();
extern gnode* graph_node;
extern unsigned num_gnodes;
extern unsigned adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned pred_adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned goal, start;
extern gnode* start_state;
extern gnode* goal_state;
extern snode* start_node;

extern unsigned long long int stat_expansions;
extern unsigned long long int stat_generated;
extern unsigned long long int minf_solution;
extern unsigned long long int stat_created;
extern unsigned long long int stat_exsarcted;
extern unsigned long long int stat_recycled;

extern unsigned solutions[MAX_SOLUTIONS][2];
extern unsigned nsolutions;
extern unsigned stat_pruned;

extern snode* snode_pool;
extern int snode_pool_idx;

extern snode* alloc_snode();
extern void initialize_parameters();
extern int backward_dijkstra(int dim);

unsigned ws_g1[MAX_WARM_STARTS];
unsigned ws_g2[MAX_WARM_STARTS];
int warm_start_count = 0;
unsigned min_f1, max_f1, min_f2, max_f2;

/* ------------------------------------------------------------------------------*/
// Generic A* for Scout Phase
// mode 1: Minimize f1 (find E1)
// mode 2: Minimize f2 (find E2)
// mode 3: Minimize Chebyshev distance (find Knee)
/* ------------------------------------------------------------------------------*/
bool run_scout_astar(int mode, unsigned *res_g1, unsigned *res_g2) {
    snode_pool_idx = 0; // Reset memory pool for this quick search
    emptyheap();
    
    // A* closed list (best scalarized f-value per node)
    double* best_f = (double*)malloc(num_gnodes * sizeof(double));
    for(int i = 0; i < num_gnodes; i++) best_f[i] = (double)LARGE;

    snode* n = alloc_snode();
    n->state = start;
    n->g1 = 0; n->g2 = 0;
    n->key = 0.0;
    insertheap(n);
    best_f[start] = 0.0;

    bool found = false;

    while (topheap() != NULL) {
        if (is_timeout()) break;
        
        n = popheap();
        if (n->state == goal) {
            *res_g1 = n->g1;
            *res_g2 = n->g2;
            found = true;
            break;
        }

        for (short d = 1; d < adjacent_table[n->state][0] * 3; d += 3) {
            unsigned nsucc = adjacent_table[n->state][d];
            unsigned cost1 = adjacent_table[n->state][d + 1];
            unsigned cost2 = adjacent_table[n->state][d + 2];
            
            unsigned newg1 = n->g1 + cost1;
            unsigned newg2 = n->g2 + cost2;
            unsigned h1 = graph_node[nsucc].h1;
            unsigned h2 = graph_node[nsucc].h2;
            
            double newkey = 0.0;
            if (mode == 1) newkey = newg1 + h1;
            else if (mode == 2) newkey = newg2 + h2;
            else if (mode == 3) {
                double n1 = (double)((newg1 + h1) - min_f1) / (max_f1 - min_f1 + 1);
                double n2 = (double)((newg2 + h2) - min_f2) / (max_f2 - min_f2 + 1);
                newkey = (n1 > n2 ? n1 : n2) * 1000000.0; // Chebyshev max
            }

            if (newkey < best_f[nsucc]) {
                best_f[nsucc] = newkey;
                snode* succ = alloc_snode();
                succ->state = nsucc;
                succ->g1 = newg1;
                succ->g2 = newg2;
                succ->key = newkey;
                insertheap(succ);
            }
        }
    }
    free(best_f);
    return found;
}

void find_warm_start_solutions() {
    warm_start_count = 0;
    min_f1 = start_state->h1;
    min_f2 = start_state->h2;
    
    unsigned res_g1, res_g2;
    
    // 1. Extreme 1 (Min f1)
    if (run_scout_astar(1, &res_g1, &res_g2)) {
        ws_g1[warm_start_count] = res_g1;
        ws_g2[warm_start_count] = res_g2;
        warm_start_count++;
        max_f2 = res_g2; // The highest f2 we should ever care about
    }
    
    // 2. Extreme 2 (Min f2)
    if (run_scout_astar(2, &res_g1, &res_g2)) {
        ws_g1[warm_start_count] = res_g1;
        ws_g2[warm_start_count] = res_g2;
        warm_start_count++;
        max_f1 = res_g1; // The highest f1 we should ever care about
    }
    
    // 3. Chebyshev Knee Point
    if (run_scout_astar(3, &res_g1, &res_g2)) {
        ws_g1[warm_start_count] = res_g1;
        ws_g2[warm_start_count] = res_g2;
        warm_start_count++;
    }
    
    // Reset pool for the actual BOA* run
    snode_pool_idx = 0;
}

/* ------------------------------------------------------------------------------*/
// Main SG-BOA* Search
/* ------------------------------------------------------------------------------*/
int sg_boastar() {
    snode** recycled_nodes = (snode**)malloc(MAX_RECYCLE * sizeof(snode*));
    int next_recycled = 0;
    nsolutions = 0;
    stat_expansions = 0;
    stat_generated = 0;
    stat_pruned = 0;

    emptyheap();
    start_node = alloc_snode();
    ++stat_created;
    start_node->state = start;
    start_node->g1 = 0; start_node->g2 = 0;
    start_node->key = 0;
    start_node->searchtree = NULL;
    insertheap(start_node);

    while (topheap() != NULL) {
        if (is_timeout()) break;

        snode* n = popheap();
        stat_exsarcted++;

        if (n->g2 >= graph_node[n->state].gmin || n->g2 + graph_node[n->state].h2 >= minf_solution) {
            stat_pruned++;
            if (next_recycled < MAX_RECYCLE) recycled_nodes[next_recycled++] = n;
            continue;
        }

        graph_node[n->state].gmin = n->g2;

        if (n->state == goal) {
            solutions[nsolutions][0] = n->g1;
            solutions[nsolutions][1] = n->g2;
            nsolutions++;
            if (minf_solution > n->g2) minf_solution = n->g2;
            continue;
        }

        ++stat_expansions;

        for (short d = 1; d < adjacent_table[n->state][0] * 3; d += 3) {
            unsigned nsucc = adjacent_table[n->state][d];
            unsigned newg1 = n->g1 + adjacent_table[n->state][d + 1];
            unsigned newg2 = n->g2 + adjacent_table[n->state][d + 2];
            unsigned h1 = graph_node[nsucc].h1;
            unsigned h2 = graph_node[nsucc].h2;

            // 1. Standard BOA* Pruning
            if (newg2 >= graph_node[nsucc].gmin || newg2 + h2 >= minf_solution) continue;

            // 2. THE MAGIC: Warm-Start Pareto Filter Pruning
            bool dominated_by_warm_start = false;
            for(int i = 0; i < warm_start_count; i++) {
                if (newg1 + h1 > ws_g1[i] && newg2 + h2 > ws_g2[i]) {
                    dominated_by_warm_start = true;
                    break;
                }
            }
            if (dominated_by_warm_start) {
                stat_pruned++;
                continue;
            }

            snode* succ;
            if (next_recycled > 0) {
                succ = recycled_nodes[--next_recycled];
                stat_recycled++;
            } else {
                succ = alloc_snode();
                ++stat_created;
            }

            succ->state = nsucc;
            stat_generated++;
            succ->searchtree = n;
            succ->g1 = newg1;
            succ->g2 = newg2;
            succ->key = (newg1 + h1) * (double)BASE + (newg2 + h2);
            insertheap(succ);
        }
    }
    free(recycled_nodes);
    return nsolutions > 0;
}

void call_sg_boastar() {
    float total_runtime, heuristic_runtime, scout_runtime, main_runtime;
    struct timeval tstart, t_heur, t_scout, tend;
    
    initialize_parameters(); // Set up start/goal pointers
    
    gettimeofday(&tstart, NULL); // Start timer

    // 1. Compute heuristics (Dijkstra)
    backward_dijkstra(1);
    backward_dijkstra(2);
    gettimeofday(&t_heur, NULL);
    
    // 2. Scout Phase (Warm-Start)
    find_warm_start_solutions();
    gettimeofday(&t_scout, NULL);
    
    // 3. Main Phase (SG-BOA*)
    sg_boastar();
    gettimeofday(&tend, NULL);

    // Calculate runtimes in seconds
    heuristic_runtime = 1.0 * (t_heur.tv_sec - tstart.tv_sec) + 1.0 * (t_heur.tv_usec - tstart.tv_usec) / 1000000.0;
    scout_runtime = 1.0 * (t_scout.tv_sec - t_heur.tv_sec) + 1.0 * (t_scout.tv_usec - t_heur.tv_usec) / 1000000.0;
    main_runtime = 1.0 * (tend.tv_sec - t_scout.tv_sec) + 1.0 * (tend.tv_usec - t_scout.tv_usec) / 1000000.0;
    total_runtime = 1.0 * (tend.tv_sec - tstart.tv_sec) + 1.0 * (tend.tv_usec - tstart.tv_usec) / 1000000.0;

    // Verify all solutions not dominated
    bool dominated = false;
    for (int i = 0; i < nsolutions; i++) {
        for (int j = 0; j < nsolutions; j++) {
            if (solutions[i][0] <= solutions[j][0] && solutions[i][1] <= solutions[j][1] && i != j) {
                dominated = true;
                break;
            }
        }
    }

    // Print Professional Summary
    printf("\n=========================================\n");
    printf("    Skeleton-Guided BOA* (SG-BOA*)       \n");
    printf("-----------------------------------------\n");
    printf("Start Node: %lld\n", start_state->id + 1);
    printf("Goal Node:  %lld\n", goal_state->id + 1);
    printf("Warm-Start Solutions Found: %d\n", warm_start_count);
    printf("Final Pareto Solutions: %d\n\n", nsolutions);

    printf("--- Execution Time (ms) ---\n");
    printf("%-25s %.3f\n", "Heuristics (Dijkstra):", heuristic_runtime * 1000);
    printf("%-25s %.3f\n", "Scout Phase (Warm-Start):", scout_runtime * 1000);
    printf("%-25s %.3f\n", "Main Search (SG-BOA*):", main_runtime * 1000);
    printf("-----------------------------------------\n");
    printf("%-25s %.3f\n\n", "Total Runtime:", total_runtime * 1000);

    printf("--- Search Statistics (Main Phase) ---\n");
    printf("%-20s %llu\n", "States Generated:", stat_generated);
    printf("%-20s %llu\n", "States Expanded:", stat_expansions);
    printf("%-20s %llu\n", "States Extracted:", stat_exsarcted);
    printf("%-20s %u\n", "States Pruned:", stat_pruned);
    printf("%-20s %llu\n", "States Created:", stat_created);
    printf("%-20s %llu\n", "States Recycled:", stat_recycled);
    printf("%-20s %s\n\n", "Contains Dominated:", dominated ? "Yes" : "No");

    printf("--- Pareto Solutions ---\n");
    for (int i = 0; i < nsolutions; i++) {
        printf("Solution %2d: g1=%-7u g2=%u\n", i + 1, solutions[i][0], solutions[i][1]);
    }
    printf("=========================================\n\n");
}