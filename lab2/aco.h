#ifndef __ACO_H__
#define __ACO_H__
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cfloat>
#include <cmath>
#include "ant.h"
using namespace std;

constexpr auto INIT_FEROMON = 0.1;

class ACO
{
	string file;							// Graph filename
	size_t problem_size = 0;				// Graph problem size
	uint32_t col_size = 10;					// Ant colony size
	double alpha = 1.0;						// Alpha param for calculating tau power in probability
	double beta = 1.0;						// Beta param for calculating eta power in probability
	double q = 1.0;							// Q param for updating pheromons
	double vapor_rate = 0.0;				// Evaporation rate param
	vector<vector<double>> graph;			// Graph path quality
	vector<vector<double>> pheromons;		// Graph path feromons
	vector<Ant> ants;						// Ants vector
	vector<uint32_t> best_path;				// Best path
	double best_path_cost = DBL_MAX;		// Best path cost
	uint32_t iterations = 10;				// Iterations of algorithm
	double a = 10.0;							// a param for tau min in MMAS algorithm

public:
	ACO(string filename) : file(filename) {}

	bool load_graph();
	void print_graph();
	void set_params(uint32_t colony_size, double alpha, double beta, double q, double a, double evaporation_rate, uint32_t iterations);
	void find_path(uint32_t start_pos, uint32_t end_pos);
	void print_best_path_and_cost();

private:
	void allocate();
	void init_pheromons();
	void realloc_init_ants();
	void evaporate();
};

#endif // __ACO_H__
