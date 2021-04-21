/*
 * Header file containing functions used in other simulations.
 */

// A Genotype is 1 ints: 0 - 5
typedef tuple<int> Genotype;
typedef Genotype State;

// A History is a sequence of States
typedef vector<State> History;

// An Ensemble is a map from Genotypes to integers (populations)
// In the notation on my blackboard and notes, the Ensemble is the set of
// occupation/population numbers for each genotype g, denoted \{n_g\}
typedef map<Genotype,uint64_t> Ensemble;
// Ensemble2 is a map from Genotypes to float numbers. Since uint64 is not large enough
typedef map<Genotype,long double> Ensemble2;

//PathEnsemble is not going to be used in the main chain
// A PathEnsemble is the equivalent for the path graph: a map from a path (a
// History) to an integer population.
//typedef map<History,uint64_t> PathEnsemble;

// Initialisation functions:
void init_genotypes(set<State> &q)
{
    for (int i=0; i<6; i++) {
        State p = make_tuple(i);
        q.emplace(p);
    }
}

void init_ends(set<State> &q)
{
    State p = make_tuple(5);
    q.emplace(p);
}

void init_ensemble(Ensemble &q, set<State> G)
{
    for (auto g = G.begin(); g != G.end(); ++g)
    {
        // for genotype in set of genotypes, set all genotype populations to zero
        q[*g] = 0;
    }
}

void reset_ensemble(Ensemble &q)
{
    for (int i=0; i<6; i++) {
        q[make_tuple(i)]=0;
    }
}

void reset_ensemble2(Ensemble2 &q)
{
    for (int i=0; i<6; i++) {
        q[make_tuple(i)]=0;
    }
}

long double outflow(State p){
    int i = get<0>(p);
    return rate4[i];
}

