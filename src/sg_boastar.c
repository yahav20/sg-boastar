/////////////////////////////////////////////////////////////////////
// Carlos Hernandez
// All rights reserved
// (Shadow-Guided BOA* implementation)
/////////////////////////////////////////////////////////////////////
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

// ---------------------------------------------------------
// SHARED GLOBALS (Defined in boastar.c and graph.c)
// ---------------------------------------------------------
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

// Shared functions from boastar.c
extern snode* alloc_snode();
extern void initialize_parameters();
extern snode* new_node();

// ---------------------------------------------------------
// SG-BOA* SPECIFIC GLOBALS
// ---------------------------------------------------------
unsigned stat_pruned_bounds = 0;
unsigned stat_pruned_lookahead = 0;

// ---------------------------------------------------------
// SG-BOA* SPECIFIC FUNCTIONS
// ---------------------------------------------------------

void sg_try_add_solution(unsigned f1, unsigned f2) {
    int insert_idx = 0;
    while (insert_idx < nsolutions && solutions[insert_idx][0] < f1) {
        insert_idx++;
    }

    if (insert_idx > 0 && solutions[insert_idx - 1][1] <= f2) {
        return; 
    }
    if (insert_idx < nsolutions && solutions[insert_idx][0] == f1 && solutions[insert_idx][1] <= f2) {
        return;
    }

    int read_idx = insert_idx;
    int write_idx = insert_idx;
    
    while (read_idx < nsolutions) {
        if (solutions[read_idx][1] >= f2) {
            read_idx++;
        } else {
            if (write_idx != read_idx) {
                solutions[write_idx][0] = solutions[read_idx][0];
                solutions[write_idx][1] = solutions[read_idx][1];
            }
            write_idx++;
            read_idx++;
        }
    }
    nsolutions = write_idx;

    if (nsolutions < MAX_SOLUTIONS) {
        for (int i = nsolutions; i > insert_idx; i--) {
            solutions[i][0] = solutions[i - 1][0];
            solutions[i][1] = solutions[i - 1][1];
        }
        solutions[insert_idx][0] = f1;
        solutions[insert_idx][1] = f2;
        nsolutions++;
    }
}

int sg_backward_dijkstra(int dim) {
    for (int i = 0; i < num_gnodes; ++i){
        graph_node[i].key = LARGE;
        graph_node[i].h1_shadow = LARGE;
        graph_node[i].h2_shadow = LARGE;
        graph_node[i].heapindex = 0;
    }
    emptyheap_dij();
    goal_state->key = 0;
    goal_state->h1_shadow = 0;
    goal_state->h2_shadow = 0;
    insertheap_dij(goal_state);

    while (topheap_dij() != NULL) {
        gnode* n = popheap_dij();
        
        if (dim == 1) n->h1 = n->key;
        else if (dim == 2) n->h2 = n->key;
        
        ++stat_expansions;
        for (short d = 1; d < pred_adjacent_table[n->id][0] * 3; d += 3) {
            gnode* pred = &graph_node[pred_adjacent_table[n->id][d]];
            unsigned cost1 = pred_adjacent_table[n->id][d + 1];
            unsigned cost2 = pred_adjacent_table[n->id][d + 2];
            
            unsigned new_key = n->key + (dim == 1 ? cost1 : cost2);
            unsigned new_shadow = (dim == 1 ? n->h2_shadow + cost2 : n->h1_shadow + cost1);

            bool improve = (pred->key > new_key) || (pred->key == new_key && (dim == 1 ? pred->h2_shadow : pred->h1_shadow) > new_shadow);

            if (improve) {
                pred->key = new_key;
                if (dim == 1) pred->h2_shadow = new_shadow;
                else pred->h1_shadow = new_shadow;
                insertheap_dij(pred);
            }
        }
    }
    return 1;
}

int sg_boastar() {
    snode** recycled_nodes = (snode**)malloc(MAX_RECYCLE * sizeof(snode*));
    if (recycled_nodes == NULL) {
        printf("FATAL: Could not allocate recycled_nodes.\n");
        exit(1);
    }
    int next_recycled = 0;
    nsolutions = 0;
    stat_pruned = 0;

    emptyheap();

    start_node = new_node();
    ++stat_created;
    start_node->state = start;
    start_node->g1 = 0;
    start_node->g2 = 0;
    start_node->key = 0;
    start_node->searchtree = NULL;
    insertheap(start_node);

    stat_expansions = 0;

    while (topheap() != NULL) {
        if (is_timeout()) {
            break; 
        }

        snode* n = popheap();
        stat_exsarcted++;
        short d;

        if (n->g2 >= graph_node[n->state].gmin || n->g2 + graph_node[n->state].h2 >= minf_solution) {
            stat_pruned++;
            if (next_recycled < MAX_RECYCLE) {
                recycled_nodes[next_recycled++] = n;
            }
            continue;
        }

        graph_node[n->state].gmin = n->g2;

        if (n->state == goal) { 
            sg_try_add_solution(n->g1, n->g2);
            if (minf_solution > n->g2) minf_solution = n->g2;
            continue;
        }

        ++stat_expansions;

        for (d = 1; d < adjacent_table[n->state][0] * 3; d += 3) {
            snode* succ;
            double newk1, newk2, newkey;
            unsigned nsucc = adjacent_table[n->state][d];
            unsigned cost1 = adjacent_table[n->state][d + 1];
            unsigned cost2 = adjacent_table[n->state][d + 2];

            unsigned newg1 = n->g1 + cost1;
            unsigned newg2 = n->g2 + cost2;
            unsigned h1 = graph_node[nsucc].h1;
            unsigned h2 = graph_node[nsucc].h2;

            if (newg2 >= graph_node[nsucc].gmin || newg2 + h2 >= minf_solution)
                continue;
            
            if (newg1 + h1 > LARGE || newg2 + h2 > LARGE) {
                stat_pruned_bounds++;
                continue;
            }

            if (graph_node[nsucc].h2_shadow != LARGE) {
                sg_try_add_solution(newg1 + h1, newg2 + graph_node[nsucc].h2_shadow);
            }
            if (graph_node[nsucc].h1_shadow != LARGE) {
                sg_try_add_solution(newg1 + graph_node[nsucc].h1_shadow, newg2 + h2);
            }

            bool dominated = false;
            if (nsolutions > 0) {
                unsigned target_f1 = newg1 + h1;
                unsigned target_f2 = newg2 + h2;
                
                int low = 0, high = nsolutions - 1;
                int best_idx = -1;
                
                while (low <= high) {
                    int mid = low + (high - low) / 2;
                    if (solutions[mid][0] <= target_f1) {
                        best_idx = mid; 
                        low = mid + 1; 
                    } else {
                        high = mid - 1;
                    }
                }
                
                if (best_idx != -1 && solutions[best_idx][1] <= target_f2) {
                    dominated = true;
                }
            }
            
            if (dominated) {
                stat_pruned_lookahead++; 
                continue; 
            }

            newk1 = newg1 + h1;
            newk2 = newg2 + h2;

            if (next_recycled > 0) { 
                succ = recycled_nodes[--next_recycled];
                stat_recycled++;
            }
            else {
                succ = new_node();
                ++stat_created;
            }

            succ->state = nsucc;
            stat_generated++;

            newkey = newk1 * (double)BASE + newk2;
            succ->searchtree = n;
            succ->g1 = newg1;
            succ->g2 = newg2;
            succ->key = newkey;
            insertheap(succ);
        }
    }

    free(recycled_nodes);
    return nsolutions > 0;
}

void call_sg_boastar()
{
    float runtime, heuristic_runtime;
    struct timeval tstart, tend;
    struct timeval compute_heuristic_time;
    initialize_parameters(); 

    gettimeofday(&tstart, NULL); 

    sg_backward_dijkstra(1);
    sg_backward_dijkstra(2);

    gettimeofday(&compute_heuristic_time, NULL); 
    
    sg_boastar();

    gettimeofday(&tend, NULL); 
    
    runtime = 1.0 * (tend.tv_sec - tstart.tv_sec) + 1.0 * (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    heuristic_runtime = 1.0 * (compute_heuristic_time.tv_sec - tstart.tv_sec) + 1.0 * (compute_heuristic_time.tv_usec - tstart.tv_usec) / 1000000.0;
    
    bool dominated = false;
    for (int i = 0; i < nsolutions; i++)
    {
        for (int j = 0; j < nsolutions; j++)
        {
            if (solutions[i][0] <= solutions[j][0] && solutions[i][1] <= solutions[j][1] && i != j)
            {
                dominated = true;
                break;
            }
        }
    }

    printf("\n=========================================\n");
    printf("          SG-BOA* Search Results           \n");
    printf("-----------------------------------------\n");
    printf("Start Node: %lld\n", start_state->id + 1);
    printf("Goal Node:  %lld\n", goal_state->id + 1);
    printf("Number of Solutions: %d\n\n", nsolutions);

    printf("--- Execution Time (ms) ---\n");
    printf("%-25s %.3f\n", "Heuristics (Dijkstra):", heuristic_runtime * 1000);
    printf("%-25s %.3f\n", "SG-BOA* Search:", (runtime - heuristic_runtime) * 1000);
    printf("-----------------------------------------\n");
    printf("%-25s %.3f\n\n", "Total Runtime:", runtime * 1000);

    printf("--- Search Statistics ---\n");
    printf("%-20s %llu\n", "States Generated:", stat_generated);
    printf("%-20s %llu\n", "States Expanded:", stat_expansions);
    printf("%-20s %llu\n", "States Extracted:", stat_exsarcted);
    printf("%-20s %u\n", "States Pruned (gmin):", stat_pruned);
    printf("%-20s %u\n", "Pruned (lookahead):", stat_pruned_lookahead);
    printf("%-20s %llu\n", "States Created:", stat_created);
    printf("%-20s %llu\n", "States Recycled:", stat_recycled);
    printf("%-20s %s\n\n", "Contains Dominated:", dominated ? "Yes" : "No");

    printf("--- Pareto Solutions ---\n");
    for (int i = 0; i < nsolutions; i++)
    {
        printf("Solution %2d: g1=%-7u g2=%u\n", i + 1, solutions[i][0], solutions[i][1]);
    }
    printf("=========================================\n\n");
}