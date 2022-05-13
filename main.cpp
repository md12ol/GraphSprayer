#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <bitset>
#include <numeric>

using namespace std;

#include "setu.h"
#include "stat.h"
#include "bitspray.h"

/*************************algorithm controls******************************/
#define NSE 30
#define alpha 0.5
#define mepl 3              //  Minimum epidemic length
#define rse 5               //  Re-try short epidemics
#define ftl ((int)NSE*10)   //  Final test length
#define verbose true
#define runs 30
#define mevs 1000
#define RIs 100
#define RE ((long)mevs/RIs)
#define popsize 50
#define verts 256
#define RNS 91207819
#define MNM 2
#define tsize 7

#define omega 0.5

#define states 12
#define Qz verts*(verts-1)/2 + 2
int Q[Qz];

bool dead[popsize];
double min_fit = 0.0;

/**************************Variable dictionary************************/
bitspray *bPop[popsize];        // Stores bitsprayers
int pop[popsize][Qz - 2];       // Stores UTAM
double fit[popsize];            // Fitness array
int b_epi[popsize];
int dx[popsize];                // Sorting index
bool mde_spread;
int fitFun;

// Variants
int var_count;
int bestEpi_varStart[max_vars];
int bestEpi_varEnd[max_vars];
vector<int> bestEpi_varProfs[max_vars];
int bestEpi_varCount;
int bestEpi_varParents[max_vars];
bitset<dna_len> bestEpi_varDNA[max_vars];
double var_prob;
int edit_lowB;
int edit_upB;
int init_bits;

// TODO: ADD THE FOLLOWING
double alpha_change = 0.05; // Try 0.05, 0.1, 0.15?
// TODO: Mild v. Severe infections (histogram of differences between variants and immunity strings)

/****************************Procedures********************************/
void initalg();                             //initialize the algorithm
void initpop();                             //initialize population
void matingevent();                         //run a mating event
void report(ostream &aus);                  //make a statistical report
void reportbest(ostream &aus, ostream &difc); //THING
void createReadMe(ostream &aus);                //tbd
void cmdLineIntro(ostream &aus);            //tbd
void cmdLineRun(int run, ostream &aus); //tbd
double fitness(int *upTri, int idx, bool final);    //tbd
double fitness(int idx, bitspray &A, bool final);   //tbd
void createUpTri(bitspray &A, int *upTri);  //tbd

/****************************Main Routine*******************************/
int main(int argc, char *argv[]) {
    /**
     * Output Root, Initial Bits, Var Prob, Edit LowB, Edit UpB, Len (0) or Spread (1)
     * nohup ./GraphSprayer "./Out/" ? ? ? ? ? &
     */

    fstream stat, best, dchar, readme, iGOut;   //statistics, best structures
    char fn[60];                                //file name construction buffer
    char *outLoc = new char[45];
    char *outRoot = argv[1];

    init_bits = (int) strtol(argv[2], nullptr, 10);
    var_prob = strtod(argv[3], nullptr);
    edit_lowB = (int) strtol(argv[4], nullptr, 10);
    edit_upB = (int) strtol(argv[5], nullptr, 10);
    //  F -> Epidemic Length; T -> Epidemic Spread
    fitFun = (int) strtol(argv[6], nullptr, 10);

    // Create directory for output
    if (fitFun == 1) {
        sprintf(outLoc, "%sOutput - ES %03dP, %02dI, %.4f%%, %02d-%02dE/",
                outRoot, popsize, init_bits, var_prob, edit_lowB, edit_upB);
    } else {
        sprintf(outLoc, "%sOutput - EL %03dP, %02dI, %.4f%%, %02d-%02dE/",
                outRoot, popsize, init_bits, var_prob, edit_lowB, edit_upB);
    }

    mkdir(outLoc, 0777);

    initalg();

    sprintf(fn, "%sbest.dat", outLoc);
    best.open(fn, ios::out);

    sprintf(fn, "%sdifc.dat", outLoc);
    dchar.open(fn, ios::out);

    sprintf(fn, "%sreadme.dat", outLoc);
    readme.open(fn, ios::out);
    createReadMe(readme);
    readme.close();

    if (verbose) {
        cmdLineIntro(cout);
    }
    if (!verbose) {
        cout << "Started" << endl;
    }

    for (int run = 0; run < runs; run++) {
        sprintf(fn, "%srun%02d.dat", outLoc, run); // File name
        stat.open(fn, ios::out);
        if (verbose)cmdLineRun(run, cout);

        initpop();

        report(stat); //report the statistics of the initial population
        for (int mev = 0; mev < mevs; mev++) {//do mating events
            matingevent();  //run a mating event
            if ((mev + 1) % RE == 0) {//Time for a report?
                if (verbose) {
                    cout << left << setw(5) << run;
                    cout << left << setw(4) << (mev + 1) / RE;
                }
                report(stat); //report statistics
            }
        }
        stat.close();

        reportbest(best, dchar);
        cout << "Done run " << run << " of " << runs - 1 << endl;
    }
    best.close();
    dchar.close();

    delete[] outLoc;
    return (0);  //keep the system happy
}

void createReadMe(ostream &aus) {
    aus << "This file contains the info about the files in this folder." << endl;
    aus << "Graph Evolution Tool." << endl;
    aus << "Fitness Function: ";
    if (mde_spread) {
        aus << "Epidemic Spread" << endl;
    } else {
        aus << "Epidemic Length" << endl;
    }
    aus << "Epidemic Model Used: SIR with Variants" << endl;
    aus << "Representation: Bitsprayers on Adjacency Matrix" << endl << endl;
    aus << "The parameter settings are as follows: " << endl;
    aus << "Number of sample epidemics: " << NSE << endl;
    aus << "Final Test Length: " << ftl << endl;
    aus << "Alpha: " << alpha << endl;
    // TODO: ADD alpha_change
    aus << "Minimum epidemic length: " << mepl << endl;
    aus << "Re-tries for short epidemics: " << rse << endl;
    aus << "Runs: " << runs << endl;
    aus << "Mating events: " << mevs << endl;
    aus << "Population size: " << popsize << endl;
    aus << "Number of vertices: " << verts << endl;
    aus << "Maximum number of mutations: " << MNM << endl;
    aus << "Tournament size: " << tsize << endl;
    aus << "Number of States: " << states << endl;
    aus << "Decay strength for diffusion characters: " << omega << endl;
    aus << "Variant Probability: " << var_prob << endl;
    aus << "DNA Length: " << dna_len << endl;
    aus << "Initial Bits Set to 1: " << init_bits << endl;
    aus << "Edit Range: " << edit_lowB << "-" << edit_upB << endl;
    aus << endl;
    aus << "The file descriptions are as follows: " << endl;
    aus << "best.dat -> the best fitness and it's associated data for each run" << endl;
    aus << "difc.dat -> the diffusion characters of the best graph for each run" << endl;
    aus << "run##.dat -> population statistics for each run" << endl;
}

void cmdLineIntro(ostream &aus) {
    aus << "Graph Evolution Tool." << endl;
    aus << "Fitness Function: ";
    if (mde_spread) {
        aus << "Epidemic Spread" << endl;
    } else {
        aus << "Epidemic Length" << endl;
    }

    aus << "Epidemic Model Used: SIR with Variants" << endl;
    aus << "Representation: Bitsprayers on Adjacency Matrix" << endl;
    aus << "Check readme.dat for more information about parameters/output." << endl;
}

void cmdLineRun(int run, ostream &aus) {
    aus << endl << "Beginning Run " << run << " of " << runs - 1 << endl;
    aus << left << setw(5) << "Run";
    aus << left << setw(4) << "RI";
    aus << left << setw(10) << "Mean";
    aus << left << setw(12) << "95% CI";
    aus << left << setw(10) << "SD";
    aus << left << setw(8) << "Best";
    aus << endl;
    aus << left << setw(5) << run;
    aus << left << setw(4) << "0";
}

void initalg() {//initialize the algorithm
    srand48(RNS);                 //read the random number seed
    for (auto &i: bPop) {
        i = new bitspray(states);
    }
}

bool necroticFilter(const int *upTri) {
    int len = Qz - 2;
    int count[2] = {0, 0};
    int bounds[2] = {1 * verts, 6 * verts};
    if (fitFun == 1) {
        bounds[0] = 1 * verts;
        bounds[1] = 2.5 * verts;
    }

    for (int i = 0; i < len; i++) {
        count[upTri[i]]++;
    }
    if (count[1] < bounds[0] || count[1] > bounds[1]) {
        return true; // DEAD
    } else {
        return false;
    }
}

double fitness(int *upTri, int idx, bool final) {//compute the fitness
    graph G(verts);      //scratch graph
    G.UTAM(upTri);
    int len, ttl;   //maximum, length, and total removed
    int cnt;             //counter for tries
    double trials[NSE];  //stores squared error for each trial
    int en;              //epidemic number
    double accu = 0.0;         //accumulator
    int best_epi = 0;
    vector<int> varBaseProf;
    var_count = 0;
    int var_parents[max_vars];
    pair<int, int> var_lens[max_vars];
    vector<int> var_profs[max_vars];
    bitset<dna_len> variants[verts];
    vector<double> v_lengths;
    v_lengths.clear();
    v_lengths.reserve(NSE);
    vector<double> v_spreads;
    v_spreads.clear();
    v_spreads.reserve(NSE);
    // TODO: ADD valphas

    if (fitFun == 0) { //  Epidemic length
        for (en = 0; en < (final ? ftl : NSE); en++) {
            cnt = 0;
            do {
                // TODO: ADD alpha_change
                G.varSIR(0, var_count, var_profs, variants,
                         var_parents, var_lens, 0.5,
                         edit_lowB, edit_upB, var_prob, init_bits);
                v_lengths.clear();
                for (int i = 0; i <= var_count; i++) {
                    v_lengths.push_back((double) (var_lens[i].second));
                }
                sort(v_lengths.begin(), v_lengths.end());
                len = (int) v_lengths.at(var_count);
                cnt++;
            } while (len < mepl && cnt < rse);
            if (final) {
                if (len > best_epi) {
                    best_epi = len;
                    for (int i = 0; i <= var_count; i++) {
                        bestEpi_varCount = var_count;
                        bestEpi_varStart[i] = var_lens[i].first;
                        bestEpi_varEnd[i] = var_lens[i].second;
                        bestEpi_varProfs[i] = var_profs[i];
                        bestEpi_varParents[i] = var_parents[i];
                        bestEpi_varDNA[i] = variants[i];
                    }
                }
                accu = -1;
            } else {
                trials[en] = len;
            }
        }
        if (!final) {
            int longest = 0;
            for (double trial: trials) {//loop over trials
                if (trial > longest) {
                    longest = (int) trial;
                }
                accu += trial;
            }
            accu = accu / NSE;
            b_epi[idx] = longest;
        }
    } else { //  Epidemic mde_spread
        for (en = 0; en < (final ? ftl : NSE); en++) {
            cnt = 0;
            do {
                // TODO: ADD alpha_change
                G.varSIR(0, var_count, var_profs, variants,
                         var_parents, var_lens, 0.5,
                         edit_lowB, edit_upB, var_prob, init_bits);
                v_lengths.clear();
                ttl = 0;
                for (int i = 0; i <= var_count; i++) {
                    v_lengths.push_back((double) (var_lens[i].second));
                    ttl += accumulate(var_profs[i].begin(), var_profs[i].end(), 0);
                }
                sort(v_lengths.begin(), v_lengths.end());
                len = (int) v_lengths.at(var_count);
                cnt++;
            } while (len < mepl && cnt < rse);
            if (final) {
                if (ttl > best_epi) {
                    best_epi = ttl;
                    for (int i = 0; i < max_vars; i++) {
                        bestEpi_varCount = var_count;
                        bestEpi_varStart[i] = var_lens[i].first;
                        bestEpi_varEnd[i] = var_lens[i].second;
                        bestEpi_varProfs[i] = var_profs[i];
                        bestEpi_varParents[i] = var_parents[i];
                        bestEpi_varDNA[i] = variants[i];
                    }
                }
                accu = -1;
            } else {
                trials[en] = ttl;
            }
        }
        if (!final) {
            int furthest = 0;
            for (double trial: trials) {//loop over trials
                if (trial > furthest) {
                    furthest = (int) trial;
                }
                accu += trial;
            }
            accu = accu / NSE;
            b_epi[idx] = furthest;
        }
    }
    return accu;  //return the fitness value
}

void createUpTri(bitspray &A, int *upTri) {//unpack the queue
    int h, t;  //head and tail of queue

    for (t = 0; t < Qz; t++) {
        Q[t] = 0;        //clear the queue
    }
    A.reset(Q, h, t);     //reset the self driving automata
    while (t < Qz - 2) {
        A.next(Q, h, t, Qz);  //run the automata
    }

    for (int i = 0; i < Qz - 2; i++) {
        upTri[i] = Q[i];
    }
}

double fitness(int idx, bitspray &A, bool final) {
    createUpTri(A, pop[idx]);
    if (necroticFilter(pop[idx])) {
        dead[idx] = true;
        return min_fit;
    }
    double fi = fitness(pop[idx], idx, final);
    return fi;
}

void printUpTri(int *upTri) {
    for (int i = 0; i < Qz - 2; i++) {
        cout << upTri[i] << " ";
    }
    cout << endl;
}

void initpop() {//initialize population
    for (int i = 0; i < popsize; i++) {
        dead[i] = false;
        do {
            bPop[i]->randomize();
            createUpTri(*bPop[i], pop[i]);
        } while (necroticFilter(pop[i]));
        fit[i] = fitness(i, *bPop[i], false);
        dx[i] = i;
    }
    min_fit = fit[0];
    for (double i: fit) {
        if (i < min_fit) {
            min_fit = i;
        }
    }
}

void matingevent() {//run a mating event
    int rp;   //loop index, random position, swap variable

    //perform tournament selection
    tselect(fit, dx, tsize, popsize); // lowest first

    bPop[dx[0]]->copy(*bPop[dx[tsize - 2]]);
    bPop[dx[1]]->copy(*bPop[dx[tsize - 1]]);
    bPop[dx[0]]->tpc(*bPop[dx[1]]);
    rp = (int) lrand48() % MNM + 1;
    bPop[dx[0]]->mutate(rp);
    rp = (int) lrand48() % MNM + 1;
    bPop[dx[1]]->mutate(rp);
    // reset dead SDAs
    dead[dx[0]] = false;
    dead[dx[1]] = false;
    fit[dx[0]] = fitness(dx[0], *bPop[dx[0]], false);
    fit[dx[1]] = fitness(dx[1], *bPop[dx[1]], false);
}

void culling() {
    for (int i = 0; i < popsize; i++) {
        if (dead[i]) {
            do {
                bPop[i]->randomize();
                createUpTri(*bPop[i], pop[i]);
            } while (necroticFilter(pop[i]));
            fit[i] = fitness(i, *bPop[i], false);
            min_fit = fit[0];
            dead[i] = false;
        }
    }
    for (double i: fit) {
        if (i < min_fit) {
            min_fit = i;
        }
    }
}

void report(ostream &aus) {//make a statistical report
    dset D;
    int deaths = 0;
    for (bool i: dead) {
        if (i) {
            deaths++;
        }
    }
    double good_fit[popsize - deaths];
    for (double &i: good_fit) {
        i = 0.0;
    }
    int cnt = 0;
    for (int i = 0; i < popsize; i++) {
        if (!dead[i]) {
            good_fit[cnt++] = fit[i];
        }
    }

    D.add(good_fit, popsize - deaths);
    int b = 0;
    for (int i = 1; i < popsize; i++) {
        if (fit[i] > fit[b]) {
            b = i; //find best fitness
        }
    }
    if (D.Rmin() != min_fit) {
        for (int i = 0; i < popsize; i++) {
            if (dead[i]) { // update min
                fit[i] = D.Rmin();
            }
        }
        min_fit = D.Rmin();
    }

    //print report
    aus << left << setw(10) << D.Rmu();
    aus << left << setw(12) << D.RCI95();
    aus << left << setw(10) << D.Rsg();
    aus << left << setw(8) << D.Rmax() << " [" << b_epi[b] << "]";
    aus << "\tDead:" << deaths << endl;
    if (verbose) {
        cout << left << setw(10) << D.Rmu();
        cout << left << setw(12) << D.RCI95();
        cout << left << setw(10) << D.Rsg();
        cout << left << setw(8) << D.Rmax() << " [" << b_epi[b] << "]";
        cout << "\tDead: " << deaths << endl;
    }
//    culling();
}

void reportbest(ostream &aus, ostream &difc) {//report the best graph
    int b;       //loop indices and best pointer
    graph G(verts);  //scratch graph
    double En;
    static double M[verts][verts];
    static double Ent[verts];

    b = 0;
    for (int i = 1; i < popsize; i++) {
        if (fit[i] > fit[b]) {
            b = i; //find best fitness
        }
    }

    //output the fitness and the associated data
    aus << fit[b] << " -fitness" << endl;

    fitness(b, *bPop[b], true);
    aus << "Epidemic Profile" << endl;
    for (int i = 0; i <= bestEpi_varCount; i++) {
        if (i == 0) {
            aus << "NA-->" << "V" << left << setw(2) << i << "\t";
        } else {
            aus << "V" << left << setw(2) << bestEpi_varParents[i] << "->V" << left << setw(2) << i << "\t";
        }
        aus << "[" << left << setw(3) << bestEpi_varStart[i] << "-";
        aus << left << setw(3) << bestEpi_varEnd[i] << "]:\t";
        for (int j = 0; j < bestEpi_varStart[i]; j++) {
            aus << "\t";
        }
        for (int j = 0; j <= bestEpi_varEnd[i] - bestEpi_varStart[i]; j++) {
            aus << bestEpi_varProfs[i].at(j) << "\t";
        }
        aus << endl;
    }
    for (int i = 0; i <= bestEpi_varCount; i++) {
        aus << "V" << i << "\t";
        aus << bestEpi_varDNA[i] << endl;
    }

    // Write the SDA
    aus << "Self-Driving Automata" << endl;
    bPop[b]->print(aus);

    G.empty(verts);
    G.UTAM(pop[b]);
    aus << "Graph" << endl;
    G.write(aus);
    aus << endl;

    for (int i = 0; i < G.size(); i++) {
        G.DiffChar(i, omega, M[i]);
    }

    for (int i = 0; i < G.size(); i++) {//loop over vertices
        En = 0.0;  //prepare En for normalizing
        for (int j = 0; j < G.size(); j++) {
            En += M[i][j];  //total the amount of gas at vertex i
        }
        for (int j = 0; j < G.size(); j++) {
            M[i][j] /= En;  //normalize so sum(gas)=1
        }
        En = 0.0;  //zero the entropy accumulator
        for (int j = 0; j < G.size();
             j++) {//build up the individual entropy terms
            if (M[i][j] > 0) {
                En += -M[i][j] * log(M[i][j]);  //this is entropy base E
            }
        }
        Ent[i] = En / log(2);  //convert entropy to Log base 2
    }

    //Now sort the entropy vector
    bool more;  //no swaps
    do {//swap out-of-order entries
        more = false;
        for (int i = 0; i < G.size() - 1; i++) {
            if (Ent[i] < Ent[i + 1]) {
                En = Ent[i];
                Ent[i] = Ent[i + 1];
                Ent[i + 1] = En;  //swap
                more = true;  //set the flag that a swap happened
            }
        }
    } while (more); //until in order

    difc << Ent[0];  //output first entropy value
    for (int i = 1; i < G.size(); i++) {
        difc << " " << Ent[i];  //output remaining values
    }
    difc << endl;  //end line
}