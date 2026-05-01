#include "include.h"
#include "boastar.h"
#include "graph.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>


void read_adjacent_table(const char* filename) {
	FILE* f;
	int i, ori, dest, dist, t;
	f = fopen(filename, "r");
	int num_arcs = 0;
	if (f == NULL) 	{
		printf("Cannot open file %s.\n", filename);
		exit(1);
	}
	fscanf(f, "%d %d", &num_gnodes, &num_arcs);
	fscanf(f, "\n");

	for (i = 0; i < num_gnodes; i++)
		adjacent_table[i][0] = 0;

	for (i = 0; i < num_arcs; i++) {
		fscanf(f, "%d %d %d %d\n", &ori, &dest, &dist, &t);
		adjacent_table[ori - 1][0]++;
		adjacent_table[ori - 1][adjacent_table[ori - 1][0] * 3 - 2] = dest - 1;
		adjacent_table[ori - 1][adjacent_table[ori - 1][0] * 3 - 1] = dist;
		adjacent_table[ori - 1][adjacent_table[ori - 1][0] * 3] = t;

		pred_adjacent_table[dest - 1][0]++;
		pred_adjacent_table[dest - 1][pred_adjacent_table[dest - 1][0] * 3 - 2] = ori - 1;
		pred_adjacent_table[dest - 1][pred_adjacent_table[dest - 1][0] * 3 - 1] = dist;
		pred_adjacent_table[dest - 1][pred_adjacent_table[dest - 1][0] * 3] = t;
	}
	fclose(f);
}

void new_graph() {
	int y;
	if (graph_node == NULL) {
		graph_node = (gnode*) calloc(num_gnodes, sizeof(gnode));
		for (y = 0; y < num_gnodes; ++y) 		{
			graph_node[y].id = y;
			graph_node[y].gmin = LARGE;
			graph_node[y].h1 = LARGE;
			graph_node[y].h2 = LARGE;
			graph_node[y].gopfirst = NULL;
		}
	}
}

// ==========================================
// Global Timeout Management
// ==========================================
double global_timeout_sec = 60.0;
struct timeval current_test_start_time = {0, 0};

bool is_timeout() {
    // If start time was never set (e.g. running standalone boa), don't timeout
    if (current_test_start_time.tv_sec == 0) {
        return false;
    }
    
    struct timeval now;
    gettimeofday(&now, NULL);
    double elapsed = (now.tv_sec - current_test_start_time.tv_sec) + 
                     (now.tv_usec - current_test_start_time.tv_usec) / 1000000.0;
    return elapsed >= global_timeout_sec;
}

void read_coordinates(const char* co_filename) {
    FILE* file = fopen(co_filename, "r");
    if (!file) {
        printf("Warning: Could not open coordinate file %s. Heuristics will be 0.\n", co_filename);
        return;
    }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        // DIMACS format: v <node_id> <lon_micro> <lat_micro>
        if (line[0] == 'v') {
            int node_id;
            int lon_micro, lat_micro;
            
            // node_id ב-DIMACS מתחיל מ-1
            if (sscanf(line, "v %d %d %d", &node_id, &lon_micro, &lat_micro) == 3) {
                // make sure to not write out of bounds if the coordinate file has extra nodes
                if (node_id - 1 < num_gnodes) {
                    graph_node[node_id - 1].lon = (double)lon_micro / 1000000.0;
                    graph_node[node_id - 1].lat = (double)lat_micro / 1000000.0;
                }
            }
        }
    }
    fclose(file);
    printf("Coordinates loaded successfully from %s\n", co_filename);
}