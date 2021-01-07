#include "aco.h"

bool ACO::load_graph()
{
	/* Method for reading graph from .aco file */

	/* Opening file */
	cout << "-------------------------------" << endl;
	cout << "Reading file: " << file << endl;
	ifstream File(file);
	if (!File.is_open())
	{
		cout << "Error while opening file" << endl;
		cout << "-------------------------------" << endl;
		return false; 
	}

	string line;
	/* Skipping comments */
	File >> line;
	while (line[0] == 'c')
	{
		getline(File, line);
		File >> line;
	}

	/* Reading problem size */
	
	if (line[0] == 'p')
	{
		File >> problem_size;
		cout << "Problem size: " << problem_size << endl;
	}

	/* Allocating memory for graph */
	allocate();

	/* Reading graph */
	string tmp;
	double value;
	for (uint32_t i = 0; i < problem_size; i++)
	{
		File >> tmp;
		for (uint32_t j = 0; j < problem_size; j++)
		{
			File >> value;
			if (value == 0)
			{
				graph[i][j] = 0;
			}
			else
			{
				graph[i][j] = 1 / value;
			}
		}
	}

	/* Closing graph file */
	File.close();

	cout << "Reading OK." << endl;
	cout << "-------------------------------" << endl;
	return true;
}

void ACO::allocate()
{
	/* Method for allocating memory for graph */
	graph.resize(problem_size);
	for (auto& g : graph)
	{
		g.resize(problem_size);
	}

}

void ACO::init_pheromons()
{
	/* Method for memory allocating and initialising for graph */
	pheromons.resize(problem_size);
	for (auto& ph : pheromons)
	{
		ph.resize(problem_size, INIT_FEROMON);
	}

}

void ACO::realloc_init_ants()
{
	/* Method for reallocing and initialising ants */
	ants.resize(col_size, Ant(problem_size, alpha, beta, q));
}

void ACO::evaporate()
{
	/* Method for evaporating pheromons */
	for (auto& phi : pheromons)
	{
		for (auto& phj : phi)
		{
			phj *= (1 - vapor_rate);
		}
	}
}

void ACO::print_graph()
{
	/* Method for printing graph */
	for (auto& i : graph)
	{
		for (auto& j : i)
		{
			cout << j << " ";
		}
		cout << endl;
	}
}

void ACO::set_params(uint32_t colony_size, double alpha, double beta, double q, double a, double evaporation_rate, uint32_t iterations)
{
	/* Method for seting additional parameters */
	col_size = colony_size;
	this->alpha = alpha;
	this->beta = beta;
	this->q = q;
	this->vapor_rate = evaporation_rate;
	this->iterations = iterations;
	this->a = a;
}

void ACO::find_path(uint32_t start_pos, uint32_t end_pos)
{
	/* Method for path finding */
	realloc_init_ants();
	init_pheromons();

	/* Constructing solution for each ant */
	double curr_path_cost;
	vector<uint32_t> best_iter_path;
	double best_iter_path_cost = DBL_MAX;
	bool best_on_iter = false;  // Flag for switching Best of all/Best on iteration (false/true)
	for (int i = 0; i < iterations; i++)
	{
		int num = 0;
		for (auto& a : ants)
		{
			//cout << "Ant #" << num << " costs: ";
			a.find_path(start_pos, end_pos, graph, pheromons);
			curr_path_cost = a.get_path_cost();
			if (curr_path_cost < best_path_cost)
			{
				best_path_cost = curr_path_cost;
				best_path = a.get_path();
			}
			if (curr_path_cost < best_iter_path_cost)
			{
				best_iter_path_cost = curr_path_cost;
				best_iter_path = a.get_path();
			}
			//cout << curr_path_cost << endl;
			if (best_on_iter)
			{
				a.UpdatePheromons(&pheromons, vapor_rate, this->a, best_iter_path_cost);
			}
			else
			{
				a.UpdatePheromons(&pheromons, vapor_rate, this->a, best_path_cost);
			}
			//a.print_path();
			a.clear();
			num++;
		}
		best_on_iter = !best_on_iter;
		evaporate();
		cout << "Iteration: " << i << endl;
		print_best_path_and_cost();
	}
}

void ACO::print_best_path_and_cost()
{
	cout << "Best path cost: " << best_path_cost << endl;
	cout << "Best path: ";
	for (auto& v : best_path)
	{
		cout << v << " ";
	}
	cout << endl;
}

