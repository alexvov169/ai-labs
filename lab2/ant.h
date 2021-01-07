#ifndef __ANT_H__
#define __ANT_H__
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

class Ant
{
	uint32_t place;					// Current ant place
	vector<uint32_t> visited;		// Visited graph vertexes
	vector<double> taueta;			// Tau power alpha multiple eta power beta ([t^a]*[n^b])
	vector<double> probabilities;	// Next travel probabilities
	double path_cost = 0;			// Summary path cost
	double alpha;					// Alpha value
	double beta;					// Beta value
	double Q;						// Q value
	uint32_t size;					// Size of area (Nki count)
public:
	Ant(uint32_t area, double alpha_val, double beta_val, double Q_val) : size(area), alpha(alpha_val), beta(beta_val), Q(Q_val) {};
	void find_path(uint32_t current_place, uint32_t end_pos, vector<vector<double>> graph, vector<vector<double>> pheromons)
	{
		place = current_place;
		visited.push_back(place);

		while (place != end_pos)
		{
			double devider = 0.0;		// Summary devider for probability finding

			/* Calculating devider and taueta vector */
			for (uint32_t i = 0; i < size; i++)
			{
				/* [t^a]*[n^b] */
				double taueta_tmp = pow(graph[place][i], alpha) * pow(pheromons[place][i], beta);
				taueta.push_back(taueta_tmp);
				devider += taueta_tmp;
			}

			/* Calculating travel probabilities vector */
			for (uint32_t i = 0; i < size; i++)
			{
				bool visited_vert = false;
				for (auto& v : visited)
				{
					if (i == v)
						visited_vert = true;
				}
				if (i != place && !visited_vert)
				{
					/* p = [t^a]*[n^b] / En([t^a]*[n^b]) */
					probabilities.push_back(taueta[i] / devider);
				}
				else
				{
					probabilities.push_back(0.0);
				}
			}

			/* Finding next travel by random from probability */
			uint32_t max_prob_place = 0;
			double max_probability = probabilities[max_prob_place];
			bool pass_ok = false;
			for (uint32_t i = 0; i < size; i++)
			{
				if (max_probability < probabilities[i])
				{
					max_prob_place = i;
					max_probability = probabilities[max_prob_place];
				}

				if (rand() % 100 < probabilities[i] * 100)
				{
					path_cost += 1 / graph[place][i];
					place = i;
					pass_ok = true;
					break;
				}
			}
			/* IF NOT OK!!! (no true statements from random) select the most probable */
			if (!pass_ok)
			{
				path_cost += 1 / graph[place][max_prob_place];
				place = max_prob_place;
			}

			/* Clearing probabilities and taueta vectors */
			taueta.clear();
			probabilities.clear();

			/* Adding current place to visited places vector */
			visited.push_back(place);
		}
	}

	void UpdatePheromons(vector<vector<double>>* pheromons, double vapor_rate, double a, double Cbest)
	{
		double delta_taubest = 1.0 / Cbest;
		double taumax = (1.0 / (1.0 - vapor_rate)) * delta_taubest;
		double taumin = taumax / a;

		for (uint32_t v = visited.size() - 1; v >= 1; v--)
		{
			(*pheromons)[visited[v]][visited[v - 1]] += delta_taubest;
			(*pheromons)[visited[v - 1]][visited[v]] += delta_taubest;

			/* If out of range [taumin, taumax] */
			if ((*pheromons)[visited[v]][visited[v - 1]] < taumin)
			{
				(*pheromons)[visited[v]][visited[v - 1]] = taumin;
			}
			if ((*pheromons)[visited[v - 1]][visited[v]] < taumin)
			{
				(*pheromons)[visited[v - 1]][visited[v]] = taumin;
			}

			if ((*pheromons)[visited[v]][visited[v - 1]] > taumax)
			{
				(*pheromons)[visited[v]][visited[v - 1]] = taumax;
			}
			if ((*pheromons)[visited[v - 1]][visited[v]] > taumax)
			{
				(*pheromons)[visited[v - 1]][visited[v]] = taumax;
			}
		}
	}

	void print_path()
	{
		for (auto& v : visited)
		{
			cout << v << " ";
		}
		cout << endl;
	}

	double get_path_cost()
	{
		return path_cost;
	}

	vector<uint32_t> get_path()
	{
		return visited;
	}

	void clear()
	{
		visited.clear();
		taueta.clear();
		probabilities.clear();
		path_cost = 0;
	}
};

#endif // __ANT_H__
