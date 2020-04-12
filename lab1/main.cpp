#include <iostream>

#include <memory>
#include <vector>
#include <tuple>

#include <algorithm>
#include <functional>
#include <cmath>

#include <random>
#include <chrono>

using namespace std;

std::default_random_engine engine(time(nullptr));

template <typename T>
class RandomInteger {
    std::uniform_int_distribution<T> distribution;
public:
    RandomInteger(T begin, T end) : distribution(begin, end) {}
    T operator()() {
        return distribution(engine);
    }
};

template <typename T>
class RandomReal {
    std::uniform_real_distribution<T> distribution;
public:
    RandomReal(T begin, T end) : distribution(begin, end) {}
    T operator()() {
        return distribution(engine);
    }
};

template <typename T>
class Random {};

template <>
class Random<int> : public RandomInteger<int> {
public:
    Random(int b, int e) : RandomInteger<int>(b, e) {}
};

template <>
class Random<size_t> : public RandomInteger<size_t> {
public:
    Random(size_t b, size_t e) : RandomInteger<size_t>(b, e) {}
};

template <>
class Random<double> : public RandomReal<double> {
public:
    Random(double b, double e) : RandomReal<double>(b, e) {}
};


template <typename T>
ostream& operator<<(ostream& stream, const unique_ptr<vector<T> >& v);

template <typename T>
ostream& operator<<(ostream& stream, const shared_ptr<vector<T> >& v);

template <typename T>
ostream& operator<<(ostream& stream, const vector<T>& v) {
    if (!v.empty()) {
        stream << "(" << v[0];
        for (auto i = v.begin()+1; i != v.end(); ++i) {
            stream << " " << *i;
        }
        stream << ")";
    }
    return stream;
}

template <typename T>
ostream& operator<<(ostream& stream, const unique_ptr<vector<T> >& v) {
    return stream << *v;
}

template <typename T>
ostream& operator<<(ostream& stream, const shared_ptr<vector<T> >& v) {
    return stream << *v;
}


typedef vector<double> fenotype;
typedef shared_ptr<fenotype> fenotype_ptr;
typedef vector<fenotype_ptr> population_type;

ostream& operator<<(ostream& stream, const vector<fenotype_ptr>& v) {
    if (!v.empty()) {
        stream << "(" << v[0];
        for (auto i = v.begin()+1; i != v.end(); ++i) {
            stream << endl << " " << *i;
        }
        stream << ")";
    }
    return stream;
}

template <typename Gen>
population_type init_population(size_t dimension_size, size_t size, Gen r) {
    population_type population;
    for (size_t i = 0; i < size; ++i) {
        fenotype f;
        for (size_t j = 0; j < dimension_size; ++j) {
            auto v = r();
            f.push_back(v);
        }
        population.push_back(fenotype_ptr(new fenotype(f)));
    }
    return population;
}

template<typename Func>
double eval_func(size_t dimension_size, Func func, const fenotype& xs);

double fitness(double y, double denominator) {
    return (1/y) / denominator;
}

bool enough_density(double delta, const population_type& p, size_t needed_close_points) {
    /**
     * sum dxi^2 < (sum dxi)^2
     * if sum dxi^2 < delta then if (sum dxi)^2 < delta then sum dxi^2 < delta
     * then we need to evaluate just (abs (sum dxi) < sqrt delta)
     **/
    /*double sqrt_delta = sqrt(delta);
    size_t needed_close_points = close_fract * p.size();
    if (needed_close_points == 0) needed_close_points = 1;
    for (size_t k = 0; k < p.size(); ++k) {
    //return any_of(p.begin(), p.end(), [](const fenotype_ptr& xs) { return xs->back() > 0.9; });
    //auto good = find_if(p.begin(), p.end(), [](const fenotype_ptr& xs) { return xs->back() > 0.9; });
    //if (good != p.end()) {
    //    size_t k = good - p.begin();
        size_t close_points_count = 0;
        for (size_t j = 0; j < k; ++j) {
            double diff = 0;
            for (size_t i = 0; i < p[k]->size(); ++i) {
                double d = (*p[k])[i] - (*p[j])[i];
                diff += d;
            }
            if (fabs(diff) < sqrt_delta) {
                ++close_points_count;
            }
        }
        for (size_t j = k; j < p.size(); ++j) {
            double diff = 0;
            for (size_t i = 0; i < p[k]->size(); ++i) {
                double d = (*p[k])[i] - (*p[j])[i];
                diff += d;
            }
            if (fabs(diff) < sqrt_delta) {
                ++close_points_count;
            }
        }
        if (close_points_count > needed_close_points) {
            return true;
        }
        //return true;
    }
    return false;*/

    if (needed_close_points == 0) needed_close_points = 1;
    for (size_t k = 0; k < p.size(); ++k) {
        size_t close_points_count = 0;
        for (size_t j = 0; j < k; ++j) {
            if (fabs(*(p[k]->end()-2) - *(p[j]->end()-2)) < delta) {
                ++close_points_count;
            }
        }
        for (size_t j = k; j < p.size(); ++j) {
            if (fabs(*(p[k]->end()-2) - *(p[j]->end()-2)) < delta) {
                ++close_points_count;
            }
        }
        if (close_points_count >= needed_close_points) {
            return true;
        }
    }
    return false;
}

population_type tournament_selection(const population_type& p, size_t tournament_size, size_t tournament_count) {
    population_type new_pop;

    if (tournament_size < 2) return p;

    for (size_t i = 0; i < tournament_count; ++i) {
        Random<size_t> random_index(0, p.size()-1);
        vector<size_t> competitors;
        for (size_t i = 0; i < tournament_size; ++i) {
            size_t index = random_index();
            if (competitors.end() == find(competitors.begin(), competitors.end(), index)) {
                competitors.push_back(index);
            } else { --i; continue; }
        }
        auto winner = competitors[0];
        for (size_t i = 1; i < competitors.size(); ++i) {
            if (p[competitors[i]]->back() > p[winner]->back()) {
                winner = competitors[i];
            }
        }
        if (new_pop.end() == find(new_pop.begin(), new_pop.end(), p[winner])) {
            new_pop.push_back(p[winner]);
        } else { --i; continue; }
    }
    return new_pop;
}

population_type linear_crossover(size_t dimension_size, fenotype xs1, fenotype xs2) {
    population_type hs({
                           fenotype_ptr(new fenotype(dimension_size)),
                           fenotype_ptr(new fenotype(dimension_size)),
                           fenotype_ptr(new fenotype(dimension_size))
                       });
    for (size_t i = 0; i < dimension_size; ++i) {
        (*hs[0])[i] = (xs1[i] + xs2[i]) / 2;
        (*hs[1])[i] = (3*xs1[i] - xs2[i]) / 2;
        (*hs[2])[i] = (-xs1[i] + 3*xs2[i]) / 2;
    }
    return hs;
}

bool fitness_greater(const fenotype_ptr& xs1, const fenotype_ptr& xs2) { return xs1->back() > xs2->back(); }
bool fitness_less(const fenotype_ptr& xs1, const fenotype_ptr& xs2) { return xs1->back() < xs2->back(); }

pair<population_type, population_type> elitist_selection(population_type& population, double elite_fract) {
    sort(population.begin(), population.end(), fitness_greater);
    size_t elites = population.size() * elite_fract;
    if (elites == 0 && population.size() != 0) elites = 1;
    return {
        population_type(population.begin(), population.begin() + elites),
        population_type(population.begin() + elites, population.end())
    };
}

population_type crossover(size_t dimension_size, const population_type& p) {
    population_type new_pop;
    Random<size_t> random_index(0, p.size()-1);

    for (size_t first = 0; first < p.size(); ++first) {
        size_t second = random_index();
        int i = 1;
        while (second == first) {
            second = random_index();
            if (i++ % 0xffffff == 0) cout << "Hanging in endless loop in crossover..." << endl;
        }
        population_type hs(linear_crossover(dimension_size, *p[first], *p[second]));
        new_pop.insert(new_pop.end(), hs.begin(), hs.end());
    }

    return new_pop;
}

population_type crossover_all(size_t dimension_size, const population_type& p, size_t max_size) {
    population_type new_pop;

    for (size_t first = 0; first < p.size() && new_pop.size() < max_size; ++first) {
        for (size_t second = 0; second < first && new_pop.size() < max_size; ++second) {
            population_type hs(linear_crossover(dimension_size, *p[first], *p[second]));
            new_pop.insert(new_pop.end(), hs.begin(), hs.end());
        }
        for (size_t second = first + 1; second < p.size() && new_pop.size() < max_size; ++second) {
            population_type hs(linear_crossover(dimension_size, *p[first], *p[second]));
            new_pop.insert(new_pop.end(), hs.begin(), hs.end());
        }
    }

    return new_pop;
}

template <typename Gen>
void mutate(size_t dimension_size, population_type& p, Gen gen) {
    for (auto& xs : p) {
        Random<size_t> random_index(0, dimension_size-1);
        (*xs)[random_index()] = gen();
    }
}

template<typename Func>
void eval_population_fitness(size_t dimension_size, Func func, population_type& population) {
    auto population_size = population.size();
    fenotype ys(population_size);
    for (size_t i = 0; i < population_size; ++i) {
        ys[i] = eval_func(dimension_size, func, *population[i]);
    }
    double denominator = 0;
    for (size_t i = 0; i < population_size; ++i) {
        denominator += 1/ys[i];
    }
    for (size_t i = 0; i < population_size; ++i) {
        population[i]->push_back(ys[i]);
        population[i]->push_back(fitness(ys[i], denominator));
    }
}

template<typename Func>
population_type genetic_algorithm(size_t dimension_size, size_t best_size, Func func, size_t population_size, double min, double max, const double delta = 1e-3) {
    Random<double> r(min, max);

    auto population = init_population(dimension_size, population_size, r);
    eval_population_fitness(dimension_size, func, population);
    const double elite_fract = .5;
    int i = 0;
    population_type bests;
    do {
//        cout << "p0 " << population.size() << "  ";
        // p
        // p/w
        // p/w * e + p/w * (1 - e)
        // p/w*e + 3 * p/w*(1 - e) = p
        // e + 3 - 3e = w
        auto winners = tournament_selection(population, 2, population.size() / (3 - 2*elite_fract));
        //cout << "g= " << (double) population.size() / (3 - 2*elite_fract) / winners.size() << endl;
        auto pair = elitist_selection(winners, elite_fract);
        auto& elite = pair.first;
//        cout << "w " << winners.size() << "  ";
        auto reproduction = crossover_all(dimension_size, elite, 3 * population.size() * (1-elite_fract) / 2);
        auto pair2 = elitist_selection(reproduction, elite_fract);
        auto& elite2 = pair2.first;
        auto& reproduction2 = pair2.second;
        mutate(dimension_size, reproduction2, r);
        population.clear();

        for (auto& i : elite) {
            i->pop_back();
            i->pop_back();
        }
        population = elite;
        population.insert(population.end(), elite2.begin(), elite2.end());
        population.insert(population.end(), reproduction2.begin(), reproduction2.end());

        eval_population_fitness(dimension_size, func, population);
        i++;
//        cout << /*population << endl <<*/ "r " << reproduction2.size() << "  e " << elite.size() << "  p " << population.size() << endl;
//        cout << *max_element(population.begin(), population.end(), fitness_less) << endl;
        if (enough_density(delta, population, population.size() - reproduction2.size()) || i > 100000) {
            bests.insert(bests.begin(), elite.begin(), elite.begin()+best_size);
            break;
        }
    } while (true);
    return bests;
}

double sphere(double x, double a, double b) {
    return a + (x+b)*(x+b);
}

template<typename Func>
void bench_genetic_algorithm(string function_name, size_t dimension_size, Func func, size_t population_size,
                             double min, double max, double delta, double expected_x, double expected_f) {
    cout << "Testing " << function_name << ":" << endl;

    const size_t best_size = 1;
    auto start = chrono::steady_clock::now();
    population_type results(genetic_algorithm(dimension_size, best_size, func, population_size, min, max, delta));
    auto end = chrono::steady_clock::now();

    for (auto& r : results) {
        r->pop_back();
    }
    cout << "Results top " << best_size << ":" << endl
         << results << endl;

    auto best = max_element(results.begin(), results.end(), fitness_less);
    cout << "Got global minimum: " << (*best)->back() << endl;
    cout << "Expected f(" << expected_x << "...) = " << expected_f << endl;

    cout << "Elapsed time: "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;

    cout << endl;
}

template<typename Func>
double eval_func(size_t dimension_size, Func func, const fenotype& xs) {
    double result = 0;
    for (size_t i = 0; i < dimension_size; ++i) {
        result += func(xs[i]);
    }
    return result;
}

class RastriginFunction {
public:
    double operator()(double x) {
        return x*x - 10*cos(2*M_PI*x);
    }
};

template<>
double eval_func(size_t d, RastriginFunction func, const fenotype& xs) {
    double result = 0;
    for (size_t i = 0; i < d; ++i) {
        result += func(xs[i]);
    }
    return 10*d + result;
}


class RosenbrockFunction {
public:
    double operator()(double x1, double x2) {
        return 1e00*(x2 - x1*x1)*(x2 - x1*x1) + (x1 - 1)*(x1 - 1);
    }
};

template<>
double eval_func(size_t dimension_size, RosenbrockFunction func, const fenotype& xs) {
    double result = 0;
    for (size_t i = 0; i < dimension_size-1; ++i) {
        result += func(xs[i], xs[i+1]);
    }
    return result;
}

class GriewankFunction {
public:
    double y1(double x) {
        return x*x/4E00;
    }
    double y2(double x, size_t i) {
        return cos(x/sqrt(i));
    }
};

template<>
double eval_func(size_t dimension_size, GriewankFunction func, const fenotype& xs) {
    double y1_acc = 0, y2_acc = 1;
    for (size_t i = 0; i < dimension_size; ++i) {
        y1_acc += func.y1(xs[i]);
        y2_acc *= func.y2(xs[i], 1+ i);
    }
    return y1_acc - y2_acc + 1;
}

struct AckleyFunction {
    double a;
    double b;
    double c;
    double y1(double x) {
        return x*x;
    }
    double y2(double x) {
        return cos(c * x);
    }
};

template<>
double eval_func(size_t d, AckleyFunction func, const fenotype& xs) {
    double y1_acc = 0, y2_acc = 0;
    for (size_t i = 0; i < d; ++i) {
        y1_acc += func.y1(xs[i]);
        y2_acc += func.y2(xs[i]);
    }
    return -func.a *
            exp(-func.b *
                sqrt((1./d) * y1_acc)) -
            exp((1./d)*y2_acc) +
            func.a + M_E;
}

int main() {

    using namespace placeholders;

    cout << "NOTE: results vectors structure: (xi..., f)" << endl;

    bench_genetic_algorithm("Sphere", 10, bind(sphere, _1, 0, 0), 300, -100, 100, 1e-320, 0, 0);
    bench_genetic_algorithm("Sphere", 10, bind(sphere, _1, 1, 2), 300, -100, 100, 1e-320, -2, 10);
    bench_genetic_algorithm("Rastrigin", 10, RastriginFunction(), 600, -5.12, 5.12, 1e-10, 0, 0);
    bench_genetic_algorithm("Rosenbrock", 10, RosenbrockFunction(), 600, -5, 10, 1e-12, 1, 0);
    bench_genetic_algorithm("Griewank", 10, GriewankFunction(), 900, -600, 600, 1e-10, 0, 0);
    bench_genetic_algorithm("Ackley", 10, AckleyFunction{20, .2, 2*M_PI}, 600, -32.768, 32.768, 1e-10, 0, 0);

    return 0;
}
