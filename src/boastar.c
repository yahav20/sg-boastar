/////////////////////////////////////////////////////////////////////
// Carlos Hernandez
// All rights reserved
/////////////////////////////////////////////////////////////////////
#include "heap.h"
#include "node.h"
#include "include.h"
#include "boastar.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdbool.h>

#define MAX_SNODES 10000000
extern bool is_timeout();

gnode* graph_node;
unsigned num_gnodes;
unsigned adjacent_table[MAXNODES][MAXNEIGH];
unsigned pred_adjacent_table[MAXNODES][MAXNEIGH];
unsigned goal, start;
gnode* start_state;
gnode* goal_state;
snode* start_node;

unsigned long long int stat_expansions = 0;
unsigned long long int stat_generated = 0;
unsigned long long int minf_solution = LARGE;
unsigned long long int stat_created = 0;
unsigned long long int stat_exsarcted = 0;
unsigned long long int stat_recycled = 0;

unsigned solutions[MAX_SOLUTIONS][2];
unsigned nsolutions = 0;
unsigned stat_pruned = 0;

snode* snode_pool = NULL;
int snode_pool_idx = 0;

snode* alloc_snode() {
    if (snode_pool == NULL) {
        snode_pool = (snode*)malloc(MAX_SNODES * sizeof(snode));
    }
    if (snode_pool_idx >= MAX_SNODES) {
        printf("FATAL: Out of snode memory!\n");
        exit(1);
    }
    return &snode_pool[snode_pool_idx++];
}

void initialize_parameters() {
    snode_pool_idx = 0;
    start_state = &graph_node[start];
    goal_state = &graph_node[goal];
    stat_percolations = 0;

    minf_solution = LARGE;
    for (int i = 0; i < num_gnodes; i++) {
        graph_node[i].gmin = LARGE;
    }
}

int backward_dijkstra(int dim) {
    for (int i = 0; i < num_gnodes; ++i){
        graph_node[i].key = LARGE;
        graph_node[i].heapindex = 0;
    }
    emptyheap_dij();
    goal_state->key = 0;
    insertheap_dij(goal_state);

    while (topheap_dij() != NULL) {
        gnode* n;
        gnode* pred;
        short d;
        n = popheap_dij();
        if (dim == 1)
            n->h1 = n->key;
        else
            n->h2 = n->key;
        ++stat_expansions;
        for (d = 1; d < pred_adjacent_table[n->id][0] * 3; d += 3) {
            pred = &graph_node[pred_adjacent_table[n->id][d]];
            int new_weight = n->key + pred_adjacent_table[n->id][d + dim];
            if (pred->key > new_weight) {
                pred->key = new_weight;
                insertheap_dij(pred);
            }
        }
    }
    return 1;
}

snode* new_node() {
    snode* state = alloc_snode();
    state->heapindex = 0;
    return state;
}

int boastar() {
    /* MAX_RECYCLE pointers on the stack (~8MB) overflow typical stack limits; use heap. */
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

        // Check for timeout at the start of each iteration to ensure we exit promptly if time is up 
        if (is_timeout()) {
            break; 
        }

        snode* n = popheap(); //best node in open
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
            //printf("GOAL [%d,%d] nsolutions:%d expanded:%llu generated:%llu heapsize:%d pruned:%d\n", n->g1, n->g2, nsolutions, stat_expansions, stat_generated, sizeheap(), stat_pruned);
            solutions[nsolutions][0] = n->g1;
            solutions[nsolutions][1] = n->g2;
            nsolutions++;
            if (nsolutions >= MAX_SOLUTIONS) {
                printf("Maximum number of solutions reached, increase MAX_SOLUTIONS!\n");
                exit(1);
            }
            if (minf_solution > n->g2)
                minf_solution = n->g2;
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

            newk1 = newg1 + h1;
            newk2 = newg2 + h2;

            if (next_recycled > 0) { //to reuse pruned nodes in memory
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

/* ------------------------------------------------------------------------------*/
// Wrapper to run BOA* and print results
/* ------------------------------------------------------------------------------*/
// Wrapper to run BOA* and print results
void call_boastar()
{
    float runtime, heuristic_runtime;
    struct timeval tstart, tend;
    unsigned long long min_cost;
    unsigned long long min_time;
    struct timeval compute_heuristic_time;
    initialize_parameters(); // Set up start/goal pointers

    gettimeofday(&tstart, NULL); // Start timer

    // Compute heuristics using backward Dijkstra
    if (backward_dijkstra(1))
        min_cost = start_state->h1;
    if (backward_dijkstra(2))
        min_time = start_state->h2;

    gettimeofday(&compute_heuristic_time, NULL); // end heuristic timer
    
    // Run BOA*
    boastar();

    gettimeofday(&tend, NULL); // End timer
    
    runtime = 1.0 * (tend.tv_sec - tstart.tv_sec) + 1.0 * (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    heuristic_runtime = 1.0 * (compute_heuristic_time.tv_sec - tstart.tv_sec) + 1.0 * (compute_heuristic_time.tv_usec - tstart.tv_usec) / 1000000.0;
    
    // Verify all solutions not dominated
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

    // Print Professional Summary
    printf("\n=========================================\n");
    printf("          BOA* Search Results            \n");
    printf("-----------------------------------------\n");
    printf("Start Node: %lld\n", start_state->id + 1);
    printf("Goal Node:  %lld\n", goal_state->id + 1);
    printf("Number of Solutions: %d\n\n", nsolutions);

    printf("--- Execution Time (ms) ---\n");
    printf("%-25s %.3f\n", "Heuristics (Dijkstra):", heuristic_runtime * 1000);
    printf("%-25s %.3f\n", "BOA* Search:", (runtime - heuristic_runtime) * 1000);
    printf("-----------------------------------------\n");
    printf("%-25s %.3f\n\n", "Total Runtime:", runtime * 1000);

    printf("--- Search Statistics ---\n");
    printf("%-20s %llu\n", "States Generated:", stat_generated);
    printf("%-20s %llu\n", "States Expanded:", stat_expansions);
    printf("%-20s %llu\n", "States Extracted:", stat_exsarcted);
    printf("%-20s %u\n", "States Pruned:", stat_pruned);
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