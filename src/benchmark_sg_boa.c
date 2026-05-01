#include "include.h"
#include "boastar.h"
#include "sg_boastar.h" 
#include "node.h"
#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>

extern gnode* start_state;
extern gnode* goal_state;

// --- Global Timeout Management ---
extern double global_timeout_sec;
extern struct timeval current_test_start_time;
extern bool is_timeout();

// --- Structs for Benchmarking ---
typedef struct {
    unsigned start;
    unsigned goal;
} Query;

typedef struct {
    bool timeout;
    int solutions;
    double time_total;
    double time_heur_ms;
    double time_scout_ms;
    double time_search_ms;
    unsigned long long generated;
    unsigned long long expanded;
    unsigned long long extracted;
    unsigned long long created;
    unsigned long long recycled;
    unsigned pruned;
} BenchStats;

// --- Time Helper ---
double get_time_ms_bench(struct timeval start, struct timeval end) {
    return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

// --- Sorting logic for Pareto solutions comparison ---
int cmp_solutions(const void* a, const void* b) {
    unsigned* solA = (unsigned*)a;
    unsigned* solB = (unsigned*)b;
    if (solA[0] != solB[0]) return solA[0] - solB[0];
    return solA[1] - solB[1];
}

/* Format definitions for the table */
#define BENCH_METRIC_W 32
#define BENCH_VAL_W    15
#define BENCH_DIFF_W   12
#define BENCH_NA       "-"   

static void print_dashes(int n) {
    for (int i = 0; i < n; i++) putchar('-');
}

void print_border(void) {
    putchar('+'); print_dashes(BENCH_METRIC_W + 2);
    putchar('+'); print_dashes(BENCH_VAL_W + 2);
    putchar('+'); print_dashes(BENCH_VAL_W + 2);
    putchar('+'); print_dashes(BENCH_DIFF_W + 2);
    putchar('+'); putchar('\n');
}

void print_row_str(const char* metric, const char* boa_val, const char* alg_val, const char* diff_val) {
    printf("| %-*s | %-*s | %-*s | %-*s |\n", 
           BENCH_METRIC_W, metric, BENCH_VAL_W, boa_val, BENCH_VAL_W, alg_val, BENCH_DIFF_W, diff_val);
}

void print_row_int_diff(const char* metric, unsigned long long boa_val, unsigned long long alg_val) {
    char boa_str[32], alg_str[32], diff_str[32];
    snprintf(boa_str, sizeof boa_str, "%llu", boa_val);
    snprintf(alg_str, sizeof alg_str, "%llu", alg_val);
    
    if (boa_val == 0 && alg_val == 0) {
        strcpy(diff_str, "0.0%");
    } else if (boa_val == 0) {
        strcpy(diff_str, "N/A");
    } else {
        double diff = ((double)alg_val - (double)boa_val) / (double)boa_val * 100.0;
        snprintf(diff_str, sizeof diff_str, "%+.1f%%", diff);
    }
    print_row_str(metric, boa_str, alg_str, diff_str);
}

void print_row_double_diff(const char* metric, double boa_val, double alg_val) {
    char boa_str[32], alg_str[32], diff_str[32];
    snprintf(boa_str, sizeof boa_str, "%.3f", boa_val);
    snprintf(alg_str, sizeof alg_str, "%.3f", alg_val);
    
    if (boa_val < 0.001) {
        strcpy(diff_str, "N/A");
    } else {
        double diff = (alg_val - boa_val) / boa_val * 100.0;
        snprintf(diff_str, sizeof diff_str, "%+.1f%%", diff);
    }
    print_row_str(metric, boa_str, alg_str, diff_str);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Usage: ./benchmark_sg_boa Maps/<map_file.txt> [max_queries=20] [timeout_secs=60]\n");
        return 1;
    }

    char* map_path = argv[1];
    int max_queries = (argc > 2) ? atoi(argv[2]) : 20; 
    global_timeout_sec = (argc > 3) ? atof(argv[3]) : 60.0; 

    // --- 1. Load the Map ---
    read_adjacent_table(map_path);
    new_graph();
    
    char co_path[256];
    strcpy(co_path, map_path);
    char* dot = strrchr(co_path, '.');
    if (dot) strcpy(dot, ".co"); else strcat(co_path, ".co");
    read_coordinates(co_path);
    
    printf("[INFO] Map loaded successfully. Nodes: %d\n", num_gnodes);

    // --- 2. Read Queries ---
    char queries_path[256];
    char prefix[256];
    strcpy(prefix, map_path);
    char* slash = strrchr(prefix, '/');
    char* map_name = slash ? slash + 1 : prefix;
    char* dash = strchr(map_name, '-');
    if (dash) *dash = '\0'; 
    sprintf(queries_path, "Queries/%s-queries", map_name);

    FILE* qfile = fopen(queries_path, "r");
    if (!qfile) {
        printf("Error: Could not open queries file %s\n", queries_path);
        return 1;
    }

    Query queries[1000];
    int num_queries = 0;
    while (fscanf(qfile, "%u %u", &queries[num_queries].start, &queries[num_queries].goal) == 2) {
        num_queries++;
        if (num_queries >= max_queries) break;
    }
    fclose(qfile);

    unsigned (*boa_sols)[2] = (unsigned (*)[2])malloc((size_t)MAX_SOLUTIONS * sizeof(*boa_sols));
    unsigned (*sg_sols)[2] = (unsigned (*)[2])malloc((size_t)MAX_SOLUTIONS * sizeof(*sg_sols));

    double total_boa_time = 0, total_sg_time = 0;
    double total_boa_search = 0, total_sg_search = 0;
    int matches = 0;

    for (int q = 0; q < num_queries; q++) {
        start = queries[q].start - 1;
        goal = queries[q].goal - 1;

        struct timeval t_start, t_heur, t_scout, t_end;

        // ==========================================
        // 1. Run Standard BOA*
        // ==========================================
        gettimeofday(&current_test_start_time, NULL); 
        gettimeofday(&t_start, NULL);
        initialize_parameters();
        backward_dijkstra(1);
        backward_dijkstra(2);
        gettimeofday(&t_heur, NULL);
        
        stat_expansions = 0; stat_generated = 0; stat_pruned = 0; 
        stat_exsarcted = 0; stat_created = 0; stat_recycled = 0; nsolutions = 0;
        boastar();
        gettimeofday(&t_end, NULL);

        BenchStats boa_stats;
        boa_stats.time_heur_ms = get_time_ms_bench(t_start, t_heur);
        boa_stats.time_search_ms = get_time_ms_bench(t_heur, t_end);
        boa_stats.time_total = get_time_ms_bench(t_start, t_end);
        boa_stats.solutions = nsolutions;
        boa_stats.generated = stat_generated;
        boa_stats.expanded = stat_expansions;
        boa_stats.pruned = stat_pruned;            
        boa_stats.extracted = stat_exsarcted;        
        boa_stats.created = stat_created;           
        boa_stats.recycled = stat_recycled;
        
        int boa_num_sols = nsolutions;
        memcpy(boa_sols, solutions, sizeof(unsigned) * 2 * (size_t)nsolutions);

        // ==========================================
        // 2. Run SG-BOA*
        // ==========================================
        gettimeofday(&current_test_start_time, NULL);
        gettimeofday(&t_start, NULL);
        initialize_parameters();
        backward_dijkstra(1);
        backward_dijkstra(2);
        gettimeofday(&t_heur, NULL);
        
        find_warm_start_solutions();
        gettimeofday(&t_scout, NULL);
        
        stat_expansions = 0; stat_generated = 0; stat_pruned = 0;
        stat_exsarcted = 0; stat_created = 0; stat_recycled = 0; nsolutions = 0;
        sg_boastar();
        gettimeofday(&t_end, NULL);

        BenchStats sg_stats;
        sg_stats.time_heur_ms = get_time_ms_bench(t_start, t_heur);
        sg_stats.time_scout_ms = get_time_ms_bench(t_heur, t_scout);
        sg_stats.time_search_ms = get_time_ms_bench(t_scout, t_end);
        sg_stats.time_total = get_time_ms_bench(t_start, t_end);
        sg_stats.solutions = nsolutions;
        sg_stats.generated = stat_generated;
        sg_stats.expanded = stat_expansions;
        sg_stats.pruned = stat_pruned;                
        sg_stats.extracted = stat_exsarcted;          
        sg_stats.created = stat_created;             
        sg_stats.recycled = stat_recycled;           

        int sg_num_sols = nsolutions;
        memcpy(sg_sols, solutions, sizeof(unsigned) * 2 * (size_t)nsolutions);

        // ==========================================
        // Verify Match
        // ==========================================
        bool is_match = (boa_num_sols == sg_num_sols);
        if (is_match) {
            qsort(boa_sols, boa_num_sols, sizeof(unsigned) * 2, cmp_solutions);
            qsort(sg_sols, sg_num_sols, sizeof(unsigned) * 2, cmp_solutions);
            for (int i = 0; i < boa_num_sols; i++) {
                if (boa_sols[i][0] != sg_sols[i][0] || boa_sols[i][1] != sg_sols[i][1]) {
                    is_match = false;
                    break;
                }
            }
        }
        if (is_match) matches++;

        // ==========================================
        // Print Table
        // ==========================================
        printf("\n--- Query %d : Start %u -> Goal %u ---\n", q + 1, start + 1, goal + 1);
        print_border();
        print_row_str("Metric", "BOA*", "SG-BOA*", "Diff %");
        print_border();
        
        print_row_str("Match Correctness?", is_match ? "YES" : "NO", is_match ? "YES" : "NO", "-");
        print_row_int_diff("Solutions Found", boa_stats.solutions, sg_stats.solutions);
        print_border();

        print_row_double_diff("Search Phase Time (ms)", boa_stats.time_search_ms, sg_stats.time_search_ms);
        print_row_double_diff("Total Runtime (ms)", boa_stats.time_total, sg_stats.time_total);
        print_border();

        print_row_int_diff("States Generated (Into Heap)", boa_stats.generated, sg_stats.generated);
        print_row_int_diff("States Expanded (Popped)", boa_stats.expanded, sg_stats.expanded);

        print_row_int_diff("States Extracted", boa_stats.extracted, sg_stats.extracted);
        print_row_int_diff("States Pruned", boa_stats.pruned, sg_stats.pruned);
        print_row_int_diff("States Created", boa_stats.created, sg_stats.created);
        print_row_int_diff("States Recycled", boa_stats.recycled, sg_stats.recycled);
        print_border();

        total_boa_time += boa_stats.time_total;
        total_sg_time += sg_stats.time_total;
        total_boa_search += boa_stats.time_search_ms;
        total_sg_search += sg_stats.time_search_ms;
    }

    // --- Final Summary ---
    printf("\n=======================================================================\n");
    printf("                  AGGREGATE SUMMARY (%d Queries)                       \n", num_queries);
    printf("=======================================================================\n");
    printf("Correctness Matches    : %d / %d (%.1f%%)\n\n", matches, num_queries, (double)matches/num_queries*100.0);
    
    printf("--- Average Times (ms) ---\n");
    printf("%-25s %-15s %-15s %-15s\n", "Metric", "BOA*", "SG-BOA*", "Improvement");
    double avg_boa_search = total_boa_search / num_queries;
    double avg_sg_search = total_sg_search / num_queries;
    printf("%-25s %-15.3f %-15.3f %+.1f%%\n", "Avg Search Time:", avg_boa_search, avg_sg_search, (avg_sg_search - avg_boa_search) / avg_boa_search * 100.0);
    
    double avg_boa_total = total_boa_time / num_queries;
    double avg_sg_total = total_sg_time / num_queries;
    printf("%-25s %-15.3f %-15.3f %+.1f%%\n", "Avg Total Runtime:", avg_boa_total, avg_sg_total, (avg_sg_total - avg_boa_total) / avg_boa_total * 100.0);
    printf("=======================================================================\n\n");

    free(boa_sols);
    free(sg_sols);
    return 0;
}