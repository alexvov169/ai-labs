#include <iostream>
#include <ctime>
#include <cmath>
#include "aco.h"

using namespace std;

class test
{
public:
	string name;				// Filename
	uint32_t start;				// Start ant place
	uint32_t end;				// End ant place
	uint32_t col_size;			// Colony size
	double alpha;				// Alpa in p = [t^a]*[n^b] / En([t^a]*[n^b])
	double beta;				// Beate in p = [t^a]*[n^b] / En([t^a]*[n^b])
	double q;					// Q in delta_taueta = Q/L
	double a;					// a param for tau min in MMAS algorithm
	double vapor_rate;			// Vaporation rate (0.0 to 1.0) where 0.0 - zero vaporating, 1.0 - full vaporating
	uint32_t iterations;		// Iterations of algorithm
	test(string n, uint32_t s, uint32_t e, uint32_t col_size, double alpha, double beta, double q, double a, double vapor_rate, uint32_t iterations) :
		name(n), 
		start(s), 
		end(e),
		col_size(col_size),
		alpha(alpha), 
		beta(beta), 
		q(q), 
		a(a),
		vapor_rate(vapor_rate),
		iterations(iterations) {}
};

/* Al - Alpha for power of 1/(connection length) */
/* Be - Beta for feromons */
/* a - param for tau min in MMAS algorithm */
/* SZ - size of colony */
/* start, end - points on graph */
/* Vapor - vaporation rate */
/* Q from [0.0, 1.0) !!!*/
/*	 "filename", start, end, SZ, Al, Be, Q, a, Vapor, iterations */
vector<test> test_files =
{
        {"ovX4.aco", 0, 8, 10, 1.0, 1.0, 0.99, 100.0, 0.1, 5},
	// {"yuzSHP55.aco", 0, 54, 10, 1.0, 1.0, 0.99, 100.0, 0.1, 5},
	// {"yuzSHP95.aco", 0, 94, 15, 1.0, 0.75, 0.99, 1000.0, 0.1, 20},
	// {"yuzSHP155.aco", 0, 154, 20, 1.0, 1.0, 0.5, 150.0, 0.1, 15}
};


int main()
{
	srand(time(NULL));
	string filename;
	cout << "Enter filename to run file or \'runall\' to run all files" << endl;
	cout << "Filename: ";
	cin >> filename;

	uint32_t start_pos;
	uint32_t end_pos;
	uint32_t col_size;
	double alpha;
	double beta;
	double q;
	double a;
	double vapor_rate;
	uint32_t iterations;
	if (filename != "runall")
	{
		cout << "Start position: ";
		cin >> start_pos;
		cout << "End position: ";
		cin >> end_pos;
		cout << "Colony size: ";
		cin >> col_size;
		cout << "Alpha: ";
		cin >> alpha;
		cout << "Beta: ";
		cin >> beta;
		cout << "Q: ";
		cin >> q;
		cout << "A param for taumin in MMAS: ";
		cin >> a;
		cout << "Vaporation rate: ";
		cin >> vapor_rate;
		cout << "Iterations count: ";
		cin >> iterations;

		test_files.clear();
		test_files.push_back(test(filename, start_pos, end_pos, col_size, alpha, beta, q, a, vapor_rate, iterations));
	}

	for (auto& i : test_files)
	{
		ACO graph(i.name);
		if (graph.load_graph())
		{
			graph.set_params(i.col_size, i.alpha, i.beta, i.q, i.a, i.vapor_rate, i.iterations);
			graph.find_path(i.start, i.end);
			cout << "-------------------------------" << endl;
			graph.print_best_path_and_cost();
		}

	}

	return 0;
}
