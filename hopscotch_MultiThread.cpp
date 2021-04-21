#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>
#include <set>
// Multithreading
#include <thread>
// Thread monitoring
#include <atomic>
#include <chrono>

/* HOPSCOTCH
 * This simulates a compound stochastic process.
 * Some beans can hop on a hopscotch court: nodes on a graph.
 * The rate of each hop is independent and tunable.
 *
 * We want to know the probability that if we have N such beans, at least one
 * bean will have reached the end of the hopscotch court (and how it did so).
 * To find this, we simulate many replicates of a game with N beans, and
 * count how many beans reached the end of the court and had different life
 * histories.
 *
 * Each state countains a number of beans. Beans can hop between different
 * states with different rates.
 * */

// Example game: court has three dimensions, each dimension is a graph.
// One dimension is a simple one-step process, the others are two-step
// processes. Each two step process has three possible states in the first
// layer, and two possible states in the end layer. Each state has possible
// "next" states in the layer above, and associated rates for these.



// User-installed headers
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


// function to return random integers from a RNG_t and a dist.
// NB: do NOT copy r

// Load parameters:
#include "params.h" // TODO load from XML file or something
using namespace std;

#include "initialisers.h"

struct simConstants {
    double ttmax;
    double dt;
    int runs;
    set<State> genotypes;
    map<State,set<State>> Neighbours;
    int mode;
};

void gillespie_beans (vector<Ensemble> *tResults, vector<Ensemble2> *tPopulations,  int seed,
        atomic<int> *progress, simConstants Params) {
    // Pull in constant parameters
    double ttmax = Params.ttmax;
    int runs = Params.runs;
    set<State> G = Params.genotypes;

    vector<long double> weights(6,0);

    // precompute weights
    for (int i=0; i<6; i++) {
        State p = make_tuple(i);
        weights[i]+= outflow(p);

        // crypt fission propensities
        if (i==2) {
                weights[i]+= rate_APClost;
        } else{
            if (i > 2) {
                weights[i] += rate_BOTHlost;
            }
        }
        //cout << "w(" << i << ") = " << weights[i] << endl;
    }

    // main simulation loop
    // Initialize RNG
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    // Perform runs simulations:
    for (int run=0; run < runs; run++) {
        // perform in a timed while loop instead? prevents stochastic hanging

		// declare dynamical variables

        double tt = 0;
        int itt = 0; // number of writes
        long double nTotal = Nbeans; // number of cells

       vector<long double> n(6,0);

        // initial condition:
        n[0]=Nbeans;
        // Sporadic cancer: all N cells initially on {0,0,0}
        long double Gamma = 0;

        vector<long double> Si(6,0);

         for (int i=0; i<6; i++) {
            Si[i] = n[i]*weights[i];
            Gamma += Si[i];
        }

        bool Continue = 1;

        int Z = runs;


        // Main loop:

        while (Continue) {
            if (tt > ttmax)
                Continue = 0;

            long double Gamma_prev = Gamma;

            long double u = gsl_ran_flat(r, 0., Gamma);

            // convert random number into event

            int ie=0;

            // determine hopping site from u

            long double remainder = u;

            while ((ie < 6) && (remainder > Si[ie])) {
                remainder -= Si[ie];
                ie++;
            }
            if (ie > 5)
                ie = 5; // correct for long double --> double precision

            remainder /= n[ie];

            // hopping
            if (remainder < rate4[ie]) {
                // make leap

                n[ie]--; // leaves this node...
                n[ie+1]++; // arrives this node

                // update propensities

                long double we = weights[ie];
                Si[ie] -= we;
                Gamma -= we;

                int it = ie + 1;
                long double wt = weights[it];
                Si[it] += wt;
                Gamma += wt;
            } else {
                // fission
                n[ie]++;
                // update propensities
                long double we = weights[ie];
                Si[ie] += we;
                Gamma += we;

                nTotal += 1.0;
            }

            // NEW: if the target is an end state (331, 341 etc), then
            // skip the rest of the simulation.

            double Deltat;
            // we want to get the population of N5. So we won't stop a run when we have at a N5 crypt.
            /*
            int IfEnd = 0;
            IfEnd += (n[5]>0);

            if (IfEnd>0) {
                Deltat = ttmax;
                tt = ttmax + 0.5;
                Continue = 0;
            } else {*/

            double Tau = 1./Gamma_prev;

            // get random time step
            Deltat = gsl_ran_exponential(r, Tau);

            // each time tt crosses 1.0, write out data
            while (tt>=(double)itt) {
                // Write out nTotal also
                (*tResults)[itt][make_tuple(-1)] = Z;

                for (int i=0; i<6; i++) {
                    if (n[i]>0) {
                        (*tResults)[itt][make_tuple(i)]+=1.0;
                        (*tPopulations)[itt][make_tuple(i)]+=n[i];
                    }
                }

                // how far through this run are we?
                double frac = tt/ttmax;
                double percent = (100*((double)run+frac)/(double)runs);
                int ipct = (int)percent;

                progress->store(ipct);

                itt++;
            }

            tt += Deltat;

            if (tt >= ttmax + 1)
                Continue = 0;
		}
    }

    // Report that we are done:
    progress->store(100);
    gsl_rng_free(r);
}



int main (int argc, char** argv) {

    // Detect number of available cores
    unsigned int nCores = thread::hardware_concurrency();
    int nThreads = nCores-4; // AGGRESSIVE RARRR
    if (nThreads < 2) nThreads = nCores;
    //cout << "nThreads" <<nThreads<< endl;

    // Get hostname id and PID
    int hostname = 0;

    // Intialize constant parameters:
    double ttmax = 80.01; // simulation time in years.
    // ^^ NB: we add a small extra amount to make sure we capture the data
    // point at exactly 80 years.
    double dt = 0.01;
    int runs ;     // runs per thread.
    int mode = 2;
    string outputfile = "";

    bool Gillespie=1;
    bool PrintKey=1;

    for (int i=1; i<argc; i++) {
        // Errors
        if ((string)argv[i] == "-o") {
            if (i == argc-1) {
                cout << "Error: no output file" << endl;
                return 1;
            } else {
                // set outputfile
                outputfile = (string)argv[i+1];
            }
        }

    if ((string)argv[i] == "--cores") {
        if (i == argc-1) {
            cout << "Error: enter number of threads" << endl;
            return 1;
        } else {
            // set number of threads
            nThreads = atoi(argv[i+1]);
        }
    }

    if ((string)argv[i] == "--runs") {
        if (i == argc-1) {
            cout << "Error: enter number of runs for each thread" << endl;
            return 1;
        } else {
            // set number of threads
            runs = atoi(argv[i+1]);
        }
    }

    if (outputfile == "") {
        cout << "Error: no output file" << endl;
        return 1;
    }


    // Flag check done

    cout << "Mode: " << mode << endl;
    cout << hostname << endl;

    // Initialize constant objects:
    // Set of states to iterate over
    set<State> G;
    init_genotypes(G);

    // Initialize output files
    ofstream output;
    output.open(outputfile);

    // Set of possible end states ([3|4][3|4]1):
    set<State> ends;
    init_ends(ends);

    // Store all possible neighbours:
    //map<State,set<State>> Neighbours = generate_neighbours(G);

    // Bundle constant parameters and objects:
    simConstants P;
    P.ttmax = ttmax;
    P.dt = dt;
    P.runs = runs;
    P.genotypes = G;
    //P.Neighbours = Neighbours;
    P.mode = mode;

    // Initialize dynamical objects:
    // Initialize results vector: vector of vector of frequencies
    // For storing simulation statistics
    int bins = (int)(ttmax+0.5+1);
    vector<vector<Ensemble>> results(nThreads);
    vector<vector<Ensemble2>>populations(nThreads);

    Ensemble tmp;
    Ensemble2 tmp2;
    reset_ensemble(tmp); // how fast is this compared to init_ensemble?
    //init_ensemble(tmp,G); // these two are equally fast
    // Add a dummy space for storing total population
    reset_ensemble2(tmp2);
    tmp[make_tuple(-1)]= 0.;
    tmp2[make_tuple(-1)]= 0.;

    for (int i=0; i<nThreads; i++) {
        for (int it=0; it<(bins+1); it++) {
            results[i].push_back(tmp);
            populations[i].push_back(tmp2);
        }
    }

    // Initialize vector of threads:
    vector<thread> vThreads(nThreads);

    // Initalize progress monitor vector:
    vector<atomic<int>> progress(nThreads);

    // Dynamical objects now initialized.
    // Run nThreads parallel simulations:
    //mode = 2
    for (int i=0; i < nThreads; i++) {
        cout << "Spawning thread " << i << "..." <<endl;
        int seed = hostname*nThreads+i;
        cout << "seed " << seed << endl;

        progress[i]=0;
        // Start simulations:
        vThreads.at(i) = thread(gillespie_beans, &results[i], &populations[i], seed,
        &progress[i], P);
    }

    // Report thread progress
    bool Report = 1;  // Set flag to true to enter loop

    time_t starttime = time(0);
    time_t endtime   = starttime+200*60*60;
    time_t now = starttime;
    time_t last_write = starttime;

    while (Report) {
        // all modes MUST be compatible with this section
        cout << "\t|";
        Report = 0; // Reset flag each loop
        double Total = 0.;

        for (auto p = progress.begin(); p != progress.end(); ++p) {
            Total += (*p).load();
            // Test whether to continue reporting:
        //    cout << (*p).load() << ", ";

            Report = (Report || (*p < 100));
        }
        //cout << "\r" << flush;
        cout << ((int)(Total/nThreads)) << "%|\r" << flush;

        this_thread::sleep_for(chrono::milliseconds(300));

        // now write to temporary file

        now = time(0);
        if (now - last_write > 60) {
            last_write = now;

            // calculate how many runs have finished
            double Total = 0;
            for (auto p = progress.begin(); p != progress.end(); ++p) {
                // the number of runs that have finished is floor(progress*nThreads/100)
                Total += floor(nThreads*((*p).load())/100);
            }
            if (Total == 0) Total = 1; // avoid division by zero at early times
            double Ztmp = runs*nThreads;//Total;

            // write out temporary results
            if (mode != 3) {
                ofstream tmpoutput;
                tmpoutput.open("timedoutput.csv");

                for (int it=0; it<bins; it++) {
                    tmpoutput << it << ", ";
                    for (auto p = G.begin(); p != G.end(); ++p) {
                        double freq = 0;

                        for (int i=0; i<nThreads; i++) {
                            freq += (double)results[i][it][*p]/Ztmp;
                        }
                        tmpoutput << freq << ", ";
                    }
                    tmpoutput << endl;
                }

                tmpoutput.close();
            }
        }
    }

    // Wait for all threads to finish
    for (int i=0; i <nThreads; i++) {
        vThreads.at(i).join();
    }

    cout << endl;
    cout << "Simulations complete." << endl;
    cout << "Writing results..." << endl;

    if (PrintKey) {
        // Write out key
        output << "tt, " ;
        for (auto p = G.begin();p != G.end(); ++p) {
            output << get<0>(*p)<< ", ";
        }
        output << endl;
    }

    // Calculate probabilities
    // With fitness, each set of results[i] has a variable number
    // of beans over time

    double Z = nThreads*runs; // we don't use any other values any more

    cout << "Accuracy: +-" << 1.0/Z << endl;

    // Write out results
    for (int it=0; it<bins; it++) {
        // probability = prob that at least one particle is on an end state
        output << it << ", ";

        double zTotal = 0;
        for (int i=0; i<nThreads; i++) {
            zTotal += results[i][it][make_tuple(-1)];
        }
        for (auto p = G.begin(); p != G.end(); ++p) {
            double freq = 0;
            for (int i=0; i < nThreads; i++) {
                // total up the results from different threads
                // Normalize by total beans in this thread at this time
                freq += (double)results[i][it][*p]/zTotal;
            }

            output << freq << ", ";
        }
        output << endl;
    }

    //write the population data
    for (int it=0; it<bins; it++) {
        output << it << ", ";

        double zTotal = 0;
        for (int i=0; i<nThreads; i++) {
            zTotal += results[i][it][make_tuple(-1)];
        }

        for (auto p = G.begin(); p != G.end(); ++p) {
            double freq = 0;
            for (int i=0; i < nThreads; i++) {
                // total up the results from different threads
                // Normalize by total beans in this thread at this time
                freq += (double)populations[i][it][*p]/zTotal;
            }

            output << freq << ", ";
        }
        output << endl;
    }



    // Close output files
    output.close();

    cout << "Finished." << endl;

    return 0;
}
