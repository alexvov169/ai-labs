#pragma once
#include <vector>
#include <iomanip>
#include <iostream>
#include <numeric>

using namespace std;

class image
{
public:
	vector<int> data;
	size_t width, heigth;
};

class hopfield_NN
{
	size_t neurons;
	vector<int> output;
	vector<vector<double>> weights;

	int calculate_neur(vector<int> input, vector<double> weight)
	{
		double value = 0.0;
		value = std::inner_product(input.begin(), input.end(), weight.begin(), 0.0);
		return value > 0.0 ? 1 : -1;
	}
public:
	hopfield_NN(size_t count)
	{
		neurons = count;

		weights.resize(count);
		for (auto& w : weights)
		{
			w.resize(count, 0.0);
		}

		output.resize(count);
	}

	vector<int> run(vector<int> input)
	{
		vector<int> next_output = input;
		int counter = 0;
		do
		{
			output = next_output;
			next_output.clear();
			int iter = 0;
			for (size_t i = 0; i < neurons; i++)
			{
				next_output.push_back(calculate_neur(output, weights[iter]));
				iter++;
			}
			counter++;
			if (counter > 1000)
				break;
		} while (next_output != output);

		return output;
	}

	void calculace_weights(vector<image> &images)
	{
		for (size_t i = 0; i < output.size(); i++)
		{
			for (size_t j = 0; j < i; j++)
			{
				double value = 0.0;
				for (auto& in : images)
				{
					value += (double)in.data[i] * (double)in.data[j];
				}
				weights[i][j] = value;
				weights[j][i] = value;
			}
		}
	}

	void print_weights()
	{
		cout << "Weights:" << endl;
		for (auto& x : weights)
		{
			for (auto& y : x)
			{
				cout << setw(2) << y << ", ";
			}
			cout << endl;
		}
	}
private:


};

