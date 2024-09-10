#include <iostream>
#include <map>
#include <string>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <chrono>
#include <thread>
#include <mutex>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/ranges>

#include "deltaGibbs.hpp"

//Disable security warnings around the use of dynamically generated format strings for Logging purposes
#pragma GCC diagnostic ignored "-Wformat-security"

//#define _TEST
//#define _TIMING

const char BuildVersion[] = "2.0.0";

//Declarations of Types

struct TargetInfo {
    std::string id;
    bool strand;
    std::string to_str() const {
        return id + ((strand) ? '-' : '+');
    }
    bool operator==(const struct TargetInfo & other){
        return this->to_str() == other.to_str();
    }
    bool operator!=(const struct TargetInfo & other){
        return this->to_str() != other.to_str();
    }
}; //strand false = fwd, true = rev

struct cmpTgtInfo {
    bool operator()(const struct TargetInfo & a, const TargetInfo & b) const {
        return a.to_str() < b.to_str();
    }
};

typedef std::set<std::string> PrimerSet;
typedef std::set<std::string> TargetSet;
typedef std::unordered_map<std::string,std::vector<struct TargetInfo>> TargetByPrimer;
typedef std::unordered_map<std::string,std::vector<std::string>> PrimerByTarget;
typedef std::map<std::string,int> FungIdxbyPrimer;
typedef std::vector<std::set<std::string *>> FungSetVec;

typedef std::unordered_map<std::string, double> BadnessMap;

enum LogLevel {ERROR = 0, WARNING = 1, INFO = 2, DEBUG = 3};

class ProgressBar{
    int max;
    int progress;
    int last;
    bool visible;
    public :
    ProgressBar(int limit = 100,bool show=true){
        max = limit;
        visible = show;
        progress = 0;
        last = 0;
    }
    void update(int step = 1){
        progress += step;
        while(visible && 100 * progress / max > last){
            std::cerr << "=";
            last++;
        }
    }
    void close(){
        if(visible){
            std::cerr << "\n";
        }
    }
};

//Declarations of Variables

const double ln2 = 0.69314718055994530941723212145818;
// See Roberts and Rosenthal (2001) Optimal Scaling for Various Met-Hastings Algorithms
// Selected as the horizontal asymptote of Fig 4 (with large proposal dimensionality)
//const double TargetAcc = 0.234;
//Based on a low dimensionallity
const double TargetAcc = 0.44;
const double AdaptAdj = 0.9;
const double InitAdapt = 1.5;
const int BATCH_COUNT = 5;
const int BATCH_SIZE = 50;
const int MIN_CONTIG_IDENT = 4;
const int MAX_AMPLICON = 1000;
const auto MaxThreads = std::thread::hardware_concurrency();
std::mutex mu;

unsigned int Seed = 0;
bool bDebug{false}, bVerbose{false}, bMemSave{false}, bCandidates(false);
int MinLen{0}, MaxLen{30}, PivotLen{50}, GenT{10000}, GenN{15000}, Threads{1};
int InitGeneration{0};
double dGMin{-12.5}, dGMax{-10.5}, Excitability{0.01};
std::string SeqFile{}, InitFile{}, OutputPrefix{"CAPDesign"}, BlacklistFile{};
std::map<char, int> NucRank{};
std::unordered_map<std::string, int> SubseqGC{};
std::unordered_set<std::string> Blacklist{};
BadnessMap SeenBadness{};
FungSetVec FungibilitySets{};
FungIdxbyPrimer p_fMap{};
std::unordered_map<std::string,long int> ExecutionTime{};


//Declaration of Functions

void init_rank_map();
void init_parser(seqan3::argument_parser & parser);
template<typename T>
std::vector<T> slicing(const std::vector<T> & v, int O, int L,bool bRC = false);
void init_primer_len();

template<class Key, class T>
bool contains(const std::map<Key,T> & m, Key key);

double calc_free_energy(std::string & primer);
std::string validate_primer(const auto & primerObj);
void load_primers(const std::string seqFile, TargetByPrimer & p_tMap, PrimerByTarget & t_pMap);
void init_subseq_gc(const TargetByPrimer & p_tMap);
std::string select_fungible_primer(PrimerSet pSet);
void run_test();
PrimerSet cover_target_set(TargetSet universe, TargetByPrimer & p_tMap, PrimerByTarget & t_pMap, const PrimerSet & pSet = PrimerSet());

std::string revcom (std::string seq);
double badness(std::string p1, std::string p2);
double evaluate_loss_function(PrimerSet pSet);
void badness_range(double * res, std::string p1, PrimerSet::iterator p2First, PrimerSet::iterator p2Last);
double badness_thread_manager(std::string p1, PrimerSet::iterator p2First, PrimerSet::iterator p2Last);
void update_loss_function(double & loss, const PrimerSet & pSet, const PrimerSet & rSet, const PrimerSet & aSet);
double p_func(double q, int gen, double temp);
template<typename ... Args>
void Log(const char * format, LogLevel lvl, Args ... args);
void init_fungible_primers(TargetByPrimer & p_tMap, PrimerByTarget & t_pMap);
void output_primer_set(PrimerSet pSet, std::ostream & os = std::cout);
void load_blacklist(std::string file);
PrimerSet load_init_set(std::string file, TargetByPrimer & p_tMap);
void output_primer_info(PrimerSet pSet,TargetByPrimer & p_tMap, PrimerByTarget &t_pMap);
void id_mandatory_primers(PrimerSet & mandatorySet,TargetSet & tSet,TargetByPrimer &p_tMap,PrimerByTarget &t_pMap);
//void mapping_blacklist(std::string file, TargetByPrimer &p_tMap,PrimerByTarget &t_pMap);


int main(int argc, char ** argv) {
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    int generation{-1}, nPInteract{0};
    double adaptiveness{InitAdapt};
    double chainHeat{1};
    double curLoss = 0;
    std::vector<double> proposalRecord{};
    //This pair of maps to vectors could be reimplemented by defining two classes
    //  Primer: with members of id and a vector of pointers to Target objects
    //  Target: with mebmers of id and strand and a vector of pointers to Primer objects
    //  Every time a primer is created it is linked to its target, and the target to its
    //  primer
    //  A std::set of Primers then is the only interface needed, rather than two global
    //  maps
    TargetByPrimer targetMap;
    PrimerByTarget primerMap;
    PrimerSet pSet{}, additionSet{}, exclusionSet{}, mandatorySet{};
    TargetSet tSet{};
    seqan3::argument_parser parser{"CAPDesign",argc,argv,seqan3::update_notifications::off};


    //Parse Command Line arguments
    init_parser(parser);
    try {
        parser.parse();
    } catch (const seqan3::argument_parser_error & ext) {
        Log("Parsing",ERROR);
        seqan3::debug_stream << ext.what() << "\n";
        throw;
    }

    //Initialize binary conversion of nucleotides
    init_rank_map();
    //Initialize potential primer lengths based on free energy parameters
    init_primer_len();
    if(bDebug){
        Log("MinLen: %d, MaxLen: %d\n",ERROR,MinLen,MaxLen);
    }
    //Handle Excitability
    if(Excitability <= 0 || Excitability > 1){
        Log("Excitability must be a non-zero value up to 1",ERROR);
        exit(1);
    }
    //Set the seed
    if(Seed == 0){
        Seed = time(NULL);
    }
    Log("Random seed set to %u",INFO,Seed);
    srand(Seed);
    if(GenT > GenN){
        Log("SA generations (%d) > than total generations (%d): All generations will be SA",WARNING,GenT,GenN);
        GenT = GenN;
    }
    if(Threads < 0){
        Log("Computation threads must be a non-negative integer: using 1 thread",WARNING);
        Threads = 1;
    } else if(Threads == 0){
        Threads = MaxThreads;
    } else if(Threads > MaxThreads){
        Log("Specified computation threads greater than system max (%d): system max",WARNING,MaxThreads);
        Threads = MaxThreads;
    }
    generation = InitGeneration-1;
    std::ofstream histStream((OutputPrefix + ".hist").c_str(),std::ios::out | std::ios::trunc);

#ifdef _TEST
    run_test();
    return 0;
#endif


    load_blacklist(BlacklistFile);
    load_primers(SeqFile,targetMap,primerMap);

    if(bCandidates){
        PrimerSet pSet{};
        for(const auto & x : targetMap){
            pSet.insert(x.first);
        }
        output_primer_set(pSet,std::cout);
        Log("Done",INFO);
        return 0;
    }

    id_mandatory_primers(mandatorySet,tSet,targetMap,primerMap);
    init_subseq_gc(targetMap);
    init_fungible_primers(targetMap,primerMap);

    //Load a user priovided initial guess, if one is provided, otherwise build one
    additionSet = load_init_set(InitFile,targetMap);
    if(additionSet.size() == 0){
        //Generate a first pass of primers in addition to the mandatory set
        additionSet = cover_target_set(tSet,targetMap,primerMap,mandatorySet);
    } else {
        //Check if the user provided set has full coverage
        //  Otherwise cover the difference
        for(const std::string & primer : additionSet){
            for(const auto & target : targetMap[primer]){
                tSet.erase(target.to_str());
            }
        }
        if(tSet.size() != 0){
            PrimerSet tmp{};
            std::set_union(mandatorySet.begin(),mandatorySet.end(),
                   additionSet.begin(), additionSet.end(),
                   std::inserter(tmp,tmp.begin()));
            tmp = cover_target_set(tSet,targetMap,primerMap,tmp);
            for(const std::string & primer : tmp){
                additionSet.insert(primer);
            }
        }
    }

    //TODO:Look into allowing subsets to be swapped out as well
    //Identify non-fungible primers and add them to the mandatory set, also remove their targets from
    //  the universe
    //This makes computation easier, and prevents the primer set size from expanding as
    //  non-fungible primers are swapped out for sets of primers covering the same targets
    for(const std::string & primer : additionSet){
        int fIdx = p_fMap[primer];
        if(mandatorySet.find(primer) != mandatorySet.end()){
            continue;
        }
        if(fIdx == -1){
            mandatorySet.insert(primer);
        } 
    }

    //Assemble the initial primer set
    std::set_union(mandatorySet.begin(),mandatorySet.end(),
                   additionSet.begin(), additionSet.end(),
                   std::inserter(pSet,pSet.begin()));
    if(pSet.size() == mandatorySet.size()){
        Log("No Primers in the minimal coverage set are swappable, Skipping simulated annealing phase",WARNING);
        output_primer_set(pSet);
        return 0;
    }

    curLoss = evaluate_loss_function(mandatorySet);
    Log("%d primers are mandatory; Minimum loss: %0.4f",INFO,mandatorySet.size(),curLoss);
    curLoss = evaluate_loss_function(pSet);

    nPInteract = pSet.size() * (pSet.size() + 1) / 2;

    chainHeat = -curLoss/2/nPInteract/log(TargetAcc);
    Log("Heating Chain...",INFO);
    ProgressBar pb(BATCH_COUNT, bVerbose || bDebug);
    int bHeating = BATCH_COUNT;
    bool bLastAdapt = false;
    while(++generation < GenN){
        double loss{curLoss}, pCrit{0};
        double rPPM = (double) (rand() % 1000000) / 1000000.0;
        PrimerSet exclusionSet{}, modSet{pSet}, fungibleSet{};
        double eCount = pSet.size() * Excitability * (1.0 - (double)(generation - 1) / (double) GenT);
        if(!bHeating && generation % 1000 == 0){
            histStream << "===Generation_" << generation << "===\n";
            output_primer_set(pSet,histStream);
        }
        std::set_difference(modSet.begin(),modSet.end(),mandatorySet.begin(),mandatorySet.end(),
                std::inserter(fungibleSet,fungibleSet.begin()));
        additionSet.clear();
        do{ //Exclude at least one primer, then keep going until eCount is hit
            std::string ePrimer = select_fungible_primer(fungibleSet);
            std::string aPrimer{};
            int r{0}, fIdx(0);
            exclusionSet.insert(ePrimer);
            modSet.erase(ePrimer);
            fungibleSet.erase(ePrimer);
            fIdx = p_fMap[ePrimer];
            do{
                r = rand() % FungibilitySets[fIdx].size();
                auto iter = FungibilitySets[fIdx].begin();
                std::advance(iter,r);
                aPrimer = *(*iter);
                //Randomly selecting from the group, until not getting the same primer 
                //E[draws] = 1/n; var = sqrt((n-1)/n^3)
                //If instead we iterated to find the index and specificially excluded that from the
                //random draw average complexity O(n/2) to find
                //In the worst case (n = 2), most of the time is at least as fast
            } while(aPrimer == ePrimer);
            additionSet.insert(aPrimer);
        } while(--eCount > 0);
        update_loss_function(loss,modSet,exclusionSet,additionSet);
        if(proposalRecord.size() >= BATCH_SIZE){
            bool increase = false;
            double mean{0};
            for(int i = 0; i < proposalRecord.size(); i++){
                mean += proposalRecord[i];
            }
            mean /= proposalRecord.size();
            if(mean < TargetAcc){
                chainHeat *= adaptiveness;
                increase = true;
            } else {
                chainHeat /= adaptiveness;
            }
            if(increase != bLastAdapt){
                bHeating--;
                pb.update();
                if(!bHeating){
                    pb.close();
                    Log("Performing Simulated Annealing...",INFO);
                }
                adaptiveness *= AdaptAdj;
                bLastAdapt = increase;
            }
            Log("BatchUpdate| m: %0.2f res: %d",DEBUG,increase);
            proposalRecord.clear();
        }
        pCrit = p_func((loss - curLoss)/nPInteract,(bHeating) ? 0 : generation,chainHeat);
        if(bHeating){
            generation--;
            if(pCrit < 1.0){
                proposalRecord.push_back(pCrit);
            }
        } else {
            Log("Gen %d | Swap: %d; o: %0.4f vs n: %0.4f; p*: %0.6f -> %s",INFO,generation+1,exclusionSet.size(),curLoss,loss,pCrit,(rPPM < pCrit) ? "Accepted" : "Rejected");
        }
        if(loss < curLoss || (generation < GenT and rPPM < pCrit)){
            curLoss = loss;
            pSet = modSet;
            for(const std::string & primer : additionSet){
                pSet.insert(primer);
            }
            
        }
    }


    std::ofstream primerStream((OutputPrefix + ".fna").c_str(),std::ios::out | std::ios::trunc);
    output_primer_set(pSet,primerStream);
    output_primer_info(pSet,targetMap,primerMap);

    Log("Done",INFO);

#ifdef _TIMING
    auto clockEnd = std::chrono::steady_clock::now();
    ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
    for(const auto & x : ExecutionTime){
        fprintf(stderr,"%25s: %15ld ns (%0.2f%%) \n",x.first.c_str(),x.second,100.0 * ((double)x.second)/((double)ExecutionTime["main"]));
    }
#endif

    return 0;
}

//Function Definitions

void init_parser(seqan3::argument_parser & parser){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    parser.add_flag(bDebug, 'D', "debug", "Include Debug Statements");
    parser.add_flag(bVerbose, 'v', "verbose", "Verbose logging");
    parser.add_flag(bMemSave, 'M', "memory-saver", "Disables recording of 'badness' of pairwise ineteractions; calculations will be repeated, but the O(n^2) structure will not be created");
    parser.add_flag(bCandidates, 'C', "output-candiates", "Outputs to stdout the complete set of candidate primers meeting energy,length, and position specifications; then exits (Useful for prefiltering primers for blacklisting purposes)");
    parser.add_option(GenT, 'a', "annealing-generations", "The number of generations until switching to gradient descent");
    parser.add_option(GenN, 'b', "total-generations", "The number of total generations for which to run");
    parser.add_option(InitGeneration, 'd', "initial-generation", "The starting generation for the annealing shedule; useful to restart or extend a run with with the --initial-primers option");
    parser.add_option(Excitability, 'e', "excitiability", "The proportion of a primer set which is swapped out at generation 0\n\t\tNote: Minimum primer swap is always 1, excitability drops as generations proceed");
    parser.add_option(dGMin, 'g', "delta-gibbs-min", "The minimum gibbs free energy for primer hybridization");
    parser.add_option(dGMax, 'G', "delta-gibbs-max", "The maximum gibbs free energy for primer hybridization");
    parser.add_option(InitFile, 'i', "initial-primers", "A fasta formatted file containing an initial primer set");
    parser.add_option(MinLen, 'l', "min-primer-len", "The min length of primer candidates to consider");
    parser.add_option(MaxLen, 'L', "max-primer-len", "The max length of primer candidates to consider");
    parser.add_option(OutputPrefix, 'o', "output-prefix", "The prefix for output files");
    parser.add_option(PivotLen, 'p', "pivot-length", "The length of target sequences");
    parser.add_option(Seed, 's', "random-seed", "A value of 0 uses local time as the seed");
    parser.add_option(Threads,'t',"threads", "The number of threads to use for calculations (0 = sys_max)");
    parser.add_option(BlacklistFile, 'x', "blacklist", "A fasta formated file containing primers to exclude from consideration");

    parser.add_positional_option(SeqFile,"A path to a fasta formated file containing targets in context");

    parser.info.author = "Zachery Dickson";
    parser.info.date = __DATE__;
    parser.info.short_description = "Program to Select a minimal set of primers for multiplex amplification of a set of target, while minimizing the number of primer dimers; the former has priority.\n Generates 4 files: *.fna, *_info_target.tab, *_info_primer.tab, and *.history; which respectively contain the primer sequences in fasta format, the primer pairs for each target, the targets and weighting for each primer, and a periodic log of sample primer sets.";
    parser.info.version = BuildVersion;
#ifdef _TIMING
    auto clockEnd = std::chrono::steady_clock::now();
    ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}


void init_rank_map(){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    int rank = 0;
    for(char nuc: {'A', 'C', 'G', 'T'}){
        NucRank[nuc] = rank;
        NucRank[tolower(nuc)] = rank;
        rank++;
    }
#ifdef _TIMING
    auto clockEnd = std::chrono::steady_clock::now();
    ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}

//Using global dGMin and Max values, calculates the primer length search space 
//First finds the trinucleotides with minimum and maximum free energy
//  Only trinucleotides with the same first and last residue need to be considered
//The max Tri will be required to be less than 0 
//MinLen = 2 * (dGMin - initGibbs) / PTriMin 
//MaxLen = 2 * (dGMax - initGibbs) / PTriMax
void init_primer_len(){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    int tmpMin{0}, tmpMax{0};
    double minFree{0}, maxFree{dGMin - 1};
    double tmpGibbs = 0;
    int size = sizeof(deltaGibbs) / sizeof(deltaGibbs[0]);
    for(int i = 0; i < size; i++){
        int j = ((i & 3) << 2) | ((i >> 2) & 3);
        double tmp = deltaGibbs[i] + deltaGibbs[j];
        if(tmp < minFree){
            minFree = tmp;
        } else if(tmp > maxFree && tmp < 0){
            if(i != 0 && i != 3 && i != 12 && i != 15){ // Contains a G or a C
                maxFree = tmp;
            }
        }
        //fprintf(stderr,"din %d,%d: %0.2f\n",i,j,tmp);
    }
    //fprintf(stderr,"min/max: %0.2f/%0.2f\n",minFree,maxFree);
    tmpMin = 2 * int((dGMin - initGibbs) / minFree) - 1;
    tmpMax = 2 * int(ceil((dGMax - initGibbs) / maxFree)) + 1;
    MinLen = (MinLen > tmpMin) ? MinLen : tmpMin;
    MaxLen = (MaxLen < tmpMax) ? MaxLen : tmpMax;
    if(MinLen < 1 || MaxLen < 1 || MinLen > MaxLen){
        Log("Could Not construct a valid candidate primer length range which meets the provided requirements for min length, max length, and allowable Gibbs free energy\n",ERROR);
    }
#ifdef _TIMING
    auto clockEnd = std::chrono::steady_clock::now();
    ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}

template<typename T>
std::vector<T> slicing(const std::vector<T> & v, int O, int L,bool bRC){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    auto first = v.begin() + O;
    auto last = v.begin() + O + L;
    std::vector<T> vector(first,last);
    if(bRC){
        std::reverse(vector.begin(),vector.end());
        for(auto & r : vector){
            r = seqan3::complement(r);
        }
    }
#ifdef _TIMING
    auto clockEnd = std::chrono::steady_clock::now();
    ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return vector;
}

template<class Key, class T>
bool contains(std::map<Key,T> & m, Key key){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    try{
        m.at(key);
    } catch (const std::out_of_range & ext){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
        return false;
    }
#ifdef _TIMING
    auto clockEnd = std::chrono::steady_clock::now();
    ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return true;
}


double calc_free_energy(std::string & primer){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    double deltaG{initGibbs};
    uint8_t dinuc{0};
    if(primer.length() < 2){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
        return deltaG;
    }
    dinuc = NucRank.at(primer[0]);
    for(int i = 1; i < primer.length();i++){
        dinuc = ((dinuc << 2) | NucRank.at(primer[i])) & 15;
        deltaG += deltaGibbs[dinuc];
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return deltaG;
}

//Takes a seqan3 'sequence' object and returns a std::string sequence
//If the primer contains invalid characters or is outside of the free energy thresholds
//an empty string is returned
std::string validate_primer(const auto & primerObj){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    std::string primerSeq{};
    double deltaG{0};
    for(int j = 0; j < primerObj.size();j++){
        if(!contains(NucRank,primerObj[j].to_char())){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
            return std::string();
        }
        primerSeq += primerObj[j].to_char();
    }
    deltaG = calc_free_energy(primerSeq);
    Log("Primer: %s ΔG: %0.4f",DEBUG,primerSeq.c_str(),deltaG);
    if(deltaG < dGMin || deltaG > dGMax || Blacklist.find(primerSeq) != Blacklist.end()){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
        return std::string();
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return primerSeq;
}


void load_primers(const std::string seqFile, TargetByPrimer & p_tMap, PrimerByTarget & t_pMap){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Loading candidate primers from %s",INFO,seqFile.c_str());
    //Open the fasta file for reading
    seqan3::sequence_file_input fin{seqFile};
    using record_type = decltype(fin)::record_type;
    //Extract Potential Primers from each target region
    for (auto & record : fin){
        int seqLen = record.sequence().size();
        //Find primers upstream of the target, then downstream
        Log("Finding primer candidates from %s",INFO,record.id().c_str());
        for(int bRC = 0; bRC < 2; bRC++){ // Process Forward then reverse primers
            int offset = (bRC) ? (seqLen + PivotLen - MinLen) / 2 : 0; 
            int sublen = (bRC) ? (seqLen - offset) : (seqLen - PivotLen + MinLen) / 2;
            struct TargetInfo target;
            target.id = record.id();
            target.id = target.id.substr(0,target.id.find(" "));
            target.strand = bRC;
            t_pMap[target.to_str()] = std::vector<std::string>();
            auto seqObj = slicing(record.sequence(),offset,sublen,bRC);
            for(int len = MinLen; len <= MaxLen; len++){ // 
                for(int i = 0; i < seqObj.size() - len + 1; i++){
                    auto primerObj = slicing(seqObj,i,len);
                    //seqan3::debug_stream << slicing(seqObj,i,len) << "\n";
                    std::string primerSeq = validate_primer(primerObj);
                    if(primerSeq.length() == 0){
                        continue;
                    }
                    try{
                        p_tMap.at(primerSeq);
                    } catch (const std::out_of_range & e){
                        p_tMap[primerSeq] = std::vector<struct TargetInfo>();
                    }
                    p_tMap[primerSeq].push_back(target);
                    t_pMap[target.to_str()].push_back(primerSeq);
                }
            }
        }
    //    seqan3::debug_stream << ">" << record.id() << "\n";
    //    seqan3::debug_stream << record.sequence() << "\n";
    //    auto r = record.sequence() | revcom;
    //    seqan3::debug_stream << r << "\n";
    //    targets.push_back(std::move(record));
    }
    //Run checks for un-targetable sequences
    std::set<std::string> toRemove{};
    for (const auto &x : t_pMap){
        if(x.second.size() == 0){
            toRemove.insert(x.first.substr(0,x.first.length()-1));
        }
    }
    if(toRemove.size() > 0){
        for(const std::string & id : toRemove){
            Log("No valid primer pairs found for %s: Removing from consideration",WARNING,id.c_str());
            struct TargetInfo target;
            target.id = id;
            target.strand = true;
            t_pMap.erase(target.to_str());
            target.strand = false;
            t_pMap.erase(target.to_str());
        }
    }
    if(bDebug){
        Log("Primer to Target Map",DEBUG);
        for (const auto & x : p_tMap){
            seqan3::debug_stream << x.first << "\n\t";
            for( const auto & info : x.second){
                seqan3::debug_stream << info.to_str() << "\t";
            }
            seqan3::debug_stream << "\n";
        }
        Log("Target to Primer Map",DEBUG);
        for (const auto & x : t_pMap){
            seqan3::debug_stream << x.first<< "\n";
            for(const auto & primer : x.second){
                seqan3::debug_stream << primer << "\t";
            } 
            seqan3::debug_stream << "\n";
        }
    }
    Log("Loaded %d primers from %d targets",INFO,p_tMap.size(),t_pMap.size() / 2);
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}

std::string revcom(std::string seq){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    int n = seq.length();
    for(int i = 0; i < n / 2; i++){
        std::swap(seq[i],seq[n - i - 1]);
    }
    for(int i = 0; i < n; i++){
        char nuc = seq[i];
        switch(nuc){
            case('A'):
                nuc = 'T';
                break;
            case('C'):
                nuc = 'G';
                break;
            case('G'):
                nuc = 'C';
                break;
            case('T'):
                nuc = 'A';
                break;
        }
        seq[i] = nuc;
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return seq;
}


//Calculate GC for all subsequences of the candidate primers
//Given the map of primers to templates, precomputes the gc content of all relevant subsequences for
//badness calculations
void init_subseq_gc(const TargetByPrimer & p_tMap){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Calculating GC for primer subsequences ...",INFO);
    ProgressBar pb(p_tMap.size(),bVerbose || bDebug);
    int progress = 0;
    int last = 0;
    for(const auto & x : p_tMap){
        pb.update();
        std::string p = x.first;
        if(p.length() < MIN_CONTIG_IDENT){
            continue;
        }
        for(int strand = 0; strand < 2; strand++,p=revcom(p)){
            for(int i = 0; i < p.length() - MIN_CONTIG_IDENT + 1; i++){
                int gc = 0;
                for(int len = 1; len < p.length() - i; len++){
                    if(p[i+len-1] == 'G' || p[i+len-1] == 'C'){
                        gc++;
                    }
                    if(len >= MIN_CONTIG_IDENT){
                        try{
                            SubseqGC[p.substr(i,len)] = gc;
                        } catch (std::out_of_range) {
                            throw;
                        }
                    }
                }
            }
        }
    }
    pb.close();
    Log("Calculated GC for primer subsequences",DEBUG);
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}

std::string select_random_primer(PrimerSet pSet, TargetByPrimer & p_tMap, PrimerByTarget & t_pMap){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Selecting a random primer",DEBUG);
    std::vector<std::string> keys{};
    int i = 0;
    //weight the chance of selecting a primer by the minimum number of primers across all targets of the primer -1 
    //target of each primer
    //i.e. 3 primers in the set {p1,p2,p3},
    //p1 -> t1, t2; p2 -> t3, t4, t5; p3 -> t6
    //t1 <- p1, p4, p5, p6; t2 <- p1; t3 <- p2,p7,p8; t4 <- p2, p7, p8, t5 <- p2, p7, p8;
    //  t6 <- p3, p9,p10
    //weight : p1 = min(|t1|,|t2|) - 1  = min(4,1) - 1 = 0,
    //         p2 = min(|t3|,|t4|,|t5|) - 1 = min(3,3,3) - 1 = 2
    //         p3 = min(|t6|) - 1 = 2 
    for(const std::string & primer : pSet){
        int weight = p_tMap.size() + 1;
        for(const struct TargetInfo & target : p_tMap[primer]){
            if(t_pMap[target.to_str()].size() < weight){
                weight = t_pMap[target.to_str()].size();
            }
        }
        if(weight < 1){
            Log("Low weight (%d) for primer %s",DEBUG,weight,primer.c_str());
        }
        keys.insert(keys.end(), weight - 1, primer);
    }
    if(keys.size() == 0){
        Log("All primers in the set are required",WARNING);
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
        return std::string("");
    }
    std::string primer = keys[rand() % keys.size()];
    Log("Selected %s",DEBUG,primer.c_str());
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return primer;
}

//Given a set of primers, and a set of primers
//  randomly selects a primer with each being weighted by the size of the fungibility set
//This weighting prevents primers in large sets from being infrequently incorporated into potential sets
//Uses the global p_fMap and FugngibilitySets variables
//It ias assumed that all primers in the set belong to a distinct, non-singular fungibility set
std::string select_fungible_primer(PrimerSet pSet){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Selecting a random fungible primer",DEBUG);
    std::vector<int> weightedIdx{};
    int idx = 0;
    std::string primer{};
    auto iter = pSet.begin();
    for(int i = 0; i < pSet.size(); i++,iter++){
        int fIdx = p_fMap[*iter];
        if(fIdx == -1){
            Log("A non-fungible primer was provided to select_fungible_primer",ERROR);
            throw "DEFINITE BUG IN THE CODE";
        }
        weightedIdx.insert(weightedIdx.end(), FungibilitySets[fIdx].size(), i);
    } 
    idx = weightedIdx[rand() % weightedIdx.size()];
    iter = pSet.begin();
    std::advance(iter,idx);
    primer = *iter;
    Log("Selected %s",DEBUG,primer.c_str());
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return primer;
}


//Greedy Set coverage: at each step of the process adds the primer which
//  maximally increases the coverage of the universe
//Returns a set of primers
//Limits search space to only primers which targets a universe restricted by the pSet 
//  The targets of pSet are expected to be absent from the universe
PrimerSet cover_target_set(TargetSet universe, TargetByPrimer & p_tMap, PrimerByTarget & t_pMap, const PrimerSet & pSet){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    PrimerSet aSet{}; //The set of primers to add to pSet;
    PrimerSet toConsider{};
    std::unordered_set<int> toIgnore{};
    ProgressBar pb(universe.size(),bVerbose || bDebug);
    std::string msg = "Selecting a minimal primer set covering %d targets";
    Log(msg.c_str(),INFO,universe.size()/2);

    //Note fungibility set indexes to ignore
    for(const auto & primer : pSet){
        int fIdx = p_fMap[primer];
        if(fIdx != -1){
            toIgnore.insert(fIdx);
        }
    }

    //Only consider one random primer per fungibility set which targets targets in the universe
    //  For non-fungible primers they themselves are considered
    for(const auto & targetStr : universe){
        for(const std::string primer : t_pMap[targetStr]){
            int fIdx = p_fMap[primer];
            if(toIgnore.find(fIdx) != toIgnore.end()){
                continue;
            }
            std::string rPrimer = primer;
            if(fIdx != -1){ //Randomly select a primer from a valid fungibility set
                int r = rand() % FungibilitySets[fIdx].size();
                auto iter = FungibilitySets[fIdx].begin();
                std::advance(iter,r);
                rPrimer = *(*iter);
            }
            toConsider.insert(rPrimer);
        }
    }

    //Main Selection Loop
    while(universe.size()){
        double bestWeight{0};
        std::vector<std::string> bestPrimerVec{};
        std::string bestPrimer{};
        PrimerSet toRemove{};
        //Find the primer(s) which maximize added coverage
        for(const std::string & primer : toConsider){
            double weight{0};
            int redundant{0};
            for(const struct TargetInfo & target : p_tMap[primer]){
                if(universe.find(target.to_str()) != universe.end()){
                    weight++;
                } else {
                    redundant++;
                }
            }
            //This ensures that redundancy is a tie breaker between otherwise equal primers
            weight = weight - 1 + 1/((double)redundant + 1); 
            if(weight >= bestWeight){
                if(weight > bestWeight){
                    bestWeight = weight;
                    bestPrimerVec.clear();
                }
                bestPrimerVec.push_back(primer);
            }
            if(weight <= 0){ //Primers adding nothing can be removed for the next iteration
                toRemove.insert(primer);
            }
        }

        if(bestWeight <= 0){
            Log("\nComplete coverage of all targets is not possible",WARNING);
            break;
        }
        pb.update(bestWeight);
        bestPrimer = bestPrimerVec[rand() % bestPrimerVec.size()];
        aSet.insert(bestPrimer);
        if(universe.size() - bestWeight == 0){
            //No need to update if the universe is covered
            break;
        }

        //Update the universe
        for(const auto & target : p_tMap[bestPrimer]){
            universe.erase(target.to_str());
        }
        //Update the consideration list
        toRemove.insert(bestPrimer);
        for(const std::string & primer : toRemove){
            toConsider.erase(primer);
        }
        Log("\n\tCurrent uncovered universe size: %d",DEBUG, universe.size()/2);
    }
    pb.close();
    Log("Found a minimal coverage set of %d primers",INFO,aSet.size() + pSet.size());
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return aSet;
}



//From Kwok et al (1990) 10.1093/nar/18.4.999
//3` AG,GA,CC -> 100 fold decrease
//3` AA -> 20 fold decrease
//if the last base is mismatched and not a T, and there is another mismatch within 4bp
//of the 3` terminus, -> 100fold decrease
double badness_oneway(std::string p1, std::string p2){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    double res{0};
    bool b3PMatch = (p1[p1.length()] == 'T');
    std::unordered_map<std::string,std::vector<int>> subseqD{};
    for(int len = MIN_CONTIG_IDENT; len < p1.length(); len++){
        for(int i = 0; i < p1.length() - len + 1; i++){
            std::string subseq = p1.substr(i, len);
            if(subseqD.find(subseq) == subseqD.end()){
                subseqD[subseq] = std::vector<int>();
            }
            subseqD[subseq].push_back(p1.length() - i - len);
        }
    }
    //Go through each subsequence of seq2 and evaluate the badness of any matches
    //Note: p2 is expected to be reverse complemented from its original state
    //  so i is actually the distance from the 3` end
    for(int len = MIN_CONTIG_IDENT; len < p2.length(); len++){
        for(int i = 0; i < p2.length() - len + 1; i++){
            std::string subseq = p2.substr(i, len);
            if(subseqD.find(subseq) != subseqD.end()){ //Found a match
                int d2 = i;
                //double base = ln2 *  (double) (len + SubseqGC[subseq]) - log((double) (d2 + 1));
                double base = pow(2.0,(double) (len + SubseqGC[subseq])) / ((double) (d2 + 1));
                for(const int & d1 : subseqD[subseq]){
                    if(!d1){
                        b3PMatch = true;
                    }
                    res +=  base / ((double)(d1 + 1));
                }
            }
        }
    }
    if(!b3PMatch){ //A non-T 3` mismatch is present
        res /= 20.0;
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return res;
}

double badness(std::string p1, std::string p2){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    double res{0};
    std::string pairStr = (p1 < p2) ? p1 + ';' + p2 : p2 + ';' + p1;
    if(!bMemSave && SeenBadness.find(pairStr) != SeenBadness.end()){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
        return SeenBadness[pairStr];
    }
    res = badness_oneway(p1,revcom(p2)) + badness_oneway(p2,revcom(p1));
    if(!bMemSave){
        mu.lock();
        SeenBadness[pairStr] = res;
        mu.unlock();
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return res;
}

void badness_range(double * res, std::string p1, PrimerSet::iterator p2First, PrimerSet::iterator p2Last){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    for(auto iter = p2First; iter != p2Last; iter++){
        *res += badness(p1,*iter);
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}

double badness_thread_manager(std::string p1, PrimerSet::iterator p2First, PrimerSet::iterator p2Last){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    double res{0};
    std::vector<std::thread *> threadList{};
    std::vector<double *> resList{};
    int remain = std::distance(p2First,p2Last);
    int blockSize = remain / Threads + 1;
    auto iterJ = p2First;
    while(remain > 0){
        std::thread * threadPtr = NULL;
        double * resPtr = new double {0};
        auto iterK = iterJ;
        if(remain < blockSize){
            iterK = p2Last;
        } else {
            std::advance(iterK,blockSize);
        }
        threadPtr = new std::thread(badness_range,resPtr,p1,iterJ,iterK);
        threadList.push_back(threadPtr);
        resList.push_back(resPtr);
        iterJ = iterK;
        remain -= blockSize;
    }
    //Synchonize
    for(std::thread *& ptr : threadList){
        ptr->join();
        delete ptr;
        ptr = NULL;
    }
    for(double *& ptr : resList){
        res += *ptr;
        delete ptr;
        ptr = NULL;
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return res;
}

double evaluate_loss_function(PrimerSet pSet){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Evaluating Loss funtion for %d primers",DEBUG,pSet.size());
    double res{0};
    for(auto iterI = pSet.begin(); iterI != pSet.end(); iterI++){
       res += badness_thread_manager(*iterI,iterI,pSet.end());
    }
    Log("Loss evaluated to be %0.4f",DEBUG,res);

#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return res;
}

//Assumes that pSet does not contain any members from either rSet or aSet
//Interactions between rSet and aSet are not removed from the loss function
void update_loss_function(double & loss, const PrimerSet & pSet, const PrimerSet & rSet, const PrimerSet & aSet){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    const PrimerSet * sets[3];
    sets[0] = &rSet;
    sets[2] = &aSet;

    //Remove the badness components from the removalSet with sign = -1
    //Add the badness components from the additionSet with sign = 1
    for(int sign = -1; sign < 2; sign +=2){
        for(auto iterI = sets[sign+1]->begin(); iterI != sets[sign+1]->end(); iterI++){
            loss += sign * badness_thread_manager(*iterI,pSet.begin(),pSet.end());
            loss += sign * badness_thread_manager(*iterI,iterI,sets[sign+1]->end());
//
//
//            //Handle interactions with primers in pSet
//            //loss += sign * badness_range(*iterI,pSet.begin(),pSet.end());
//            for(auto iterJ = pSet.begin(); iterJ != pSet.end(); iterJ++){
//                loss += sign * badness(*iterI,*iterJ);
//            }
//            //loss += sign * badness_range(*iterI,iterI,sets[sign+1]->end());
//            //Handle interactions within the a/rSet
//            for(auto iterK = iterI; iterK != sets[sign+1]->end();iterK++){
//                loss += sign * badness(*iterI,*iterK);
//            } 
        }
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}


double p_func(double q,int gen,double temp){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    if(q < 0){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
        return 1.0;
    }
    if(gen >= GenT){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
        return 0.0;
    }
    double res (0);
    Log("Calculating p_func for ΔLoss exp(-%0.4f * (1+log_10(%d)) / %0.4f)",DEBUG,q,gen+1,temp);
    res = exp(-q *  (1+log10((double)(gen+1))) / temp);
    Log("Calculated p to be %0.4f", DEBUG,res);
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return res;
}



//Determines which primers have the same set of targets
//We can proceed target by target and test all primers for the target rather than comparing all primers
//Σn_i^2 vs N^2 where all n are << N;
void init_fungible_primers(TargetByPrimer & p_tMap, PrimerByTarget & t_pMap){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Determining primer fungibility...",INFO);
    ProgressBar pb(t_pMap.size(),bDebug || bVerbose);
    for(auto & x : t_pMap){
        pb.update();
        int idx = FungibilitySets.size();
        for(std::string & primer : x.second){
            //skip primers already assigned to fungability groups
            if(p_fMap.find(primer) != p_fMap.end()){
                continue;
            }
            //As processing is done on a per target basis and fungibility is defined by targets
            //  Only new fungibility groups created for this target need to be examined
            //As all primers withing a fungability set are interchangable, a comparison need only be
            //  made to the first primer in each set
            bool bFung = 0;
            for(int i = idx; i < FungibilitySets.size(); i++){
                bFung = 1;
                std::string * primerF = *(FungibilitySets[i].begin());
                if(p_tMap[primer].size() != p_tMap[*primerF].size()){
                    continue;
                }
                for(int k = 0; k < p_tMap[primer].size() && bFung; k++){
                    if(p_tMap[primer][k] != p_tMap[*primerF][k]){
                        bFung = 0;
                    }
                }
                if(bFung){
                    FungibilitySets[i].insert(&primer);
                    p_fMap[primer] = -1;
                    break;
                }
            }
            //This will occur if the loop above is never entered (first novel primer for the target),
            //  or the primer is not fungible with any other primer for the target so far
            if(!bFung){
                p_fMap[primer] = -1;
                FungibilitySets.push_back(std::set<std::string *>{&primer});
            }
        }
        //Remove singular sets (i.e. non-fungible)
        // and index the fungible sets
        auto iter =  FungibilitySets.begin();
        //auto iterE = FungibilitySets.end();
        std::advance(iter,idx);
        FungibilitySets.erase(
        std::remove_if(iter,FungibilitySets.end(),
                    [](const std::set<std::string *> & x){return x.size() < 2;}) /*;//*/,
                FungibilitySets.end());
    }
    for(int i = 0; i < FungibilitySets.size(); i++){
        for(std::string * primer : FungibilitySets[i]){
            p_fMap[*primer] = i;
        }
    }
    FungibilitySets.shrink_to_fit();
    pb.close();
    Log("Initialized Primer Fungibility",DEBUG);
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}


void output_primer_set(PrimerSet pSet, std::ostream & os){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Logging primers ...",INFO);
    int counter = 0;
    for(const auto & primer : pSet){
        char buf[100];
        sprintf(buf,">Primer_%04d",counter++);
        os << buf << "\n" << primer << "\n";
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}


void load_blacklist(std::string file){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    if(file.empty()){
        return;
    }
    Log("Loading blacklist from %s ...",INFO,file.c_str());
    seqan3::sequence_file_input fin{file};
    using record_type = decltype(fin)::record_type;
    //Extract Potential Primers from each target region
    for (auto & record : fin){
        std::string primerSeq{};
        for(int j = 0; j < record.sequence().size();j++){
            primerSeq += record.sequence()[j].to_char();
        }
        Blacklist.insert(primerSeq);
    }
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}

PrimerSet load_init_set(std::string file, TargetByPrimer & p_tMap){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    PrimerSet pSet{};
    if(file.empty()){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
        return pSet;
    }
    Log("Loading initial primer set from %s ...",INFO,file.c_str());
    seqan3::sequence_file_input fin{file};
    using record_type = decltype(fin)::record_type;
    //Extract Potential Primers from each target region
    for (auto & record : fin){
        std::string primerSeq{};
        for(int j = 0; j < record.sequence().size();j++){
            if(!contains(NucRank,record.sequence()[j].to_char())){
                Log("Skipping primer containing non-standard nucleotides",WARNING);
                goto skipSeq;
            }
            primerSeq += record.sequence()[j].to_char();
        }
        if(p_tMap.find(primerSeq) == p_tMap.end()){
            Log("Skipping invalid primer (%s)",WARNING,primerSeq.c_str());
            continue;
        }
        pSet.insert(primerSeq);
        skipSeq:;
    }
    Log("Loaded initial set of %d primers",DEBUG,pSet.size());
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
    return pSet;
}


void output_primer_info(PrimerSet pSet,TargetByPrimer & p_tMap, PrimerByTarget &t_pMap){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Logging primer and target metadata ...",INFO);
    std::ofstream primerInfo((OutputPrefix + "_info_primer.tsv").c_str(),std::ios::out | std::ios::trunc);
    std::ofstream targetInfo((OutputPrefix + "_info_target.tsv").c_str(),std::ios::out | std::ios::trunc);
    int nextID = 0;
    primerInfo << "Primer\tScore\tWeight\tFwdTargets\tRevTargets\n";
    for(const std::string & primer : pSet){
        char id[100], score[1000];
        sprintf(id,"Primer_%04d",nextID++);
        std::string strandList[2];
        double weight{0.0};
        double bd{0.0};
        for(const std::string & other : pSet){
            bd += badness(primer,other);
        }
        for(const struct TargetInfo & target : p_tMap.at(primer)){
            if(!strandList[target.strand].empty()){
                strandList[target.strand] += ",";
            }
            strandList[target.strand] += target.id;
            int counter = 0;
            for(const std::string & subprimer : t_pMap.at(target.to_str())){
                if(pSet.find(subprimer) != pSet.end()){
                    counter++;
                }
            }
            if(counter == 0){
                Log("Primer %s not found in primer set for target %s",ERROR,primer.c_str(),target.to_str());
                throw "DEFINITE BUG";
            } else {
                weight += 1.0 / (double)counter;
            }
        }
        primerInfo << id << "\t" << (int)bd << "\t"<< weight << "\t" << strandList[0] << "\t" << strandList[1] << "\n";
    }
    primerInfo.close();
    targetInfo << "Target\tPairs\tFwdPrimers\tRevPrimers\n";
    auto iter0 = pSet.begin();
    for(auto & x : t_pMap){
        std::string id = x.first.substr(0,x.first.length()-1);
        std::string strand = x.first.substr(x.first.length()-1,1);
        int nPrimer[] = {0,0};
        if(strand == "-"){ // All handling will be done based on the forward target
            continue;
        }
        std::vector<std::string> * sets[2];
        sets[0] = &x.second;
        sets[1] = &t_pMap[id + "-"];
        std::string strandList[2];
        for(int s = 0; s < 2; s++){ 
            for(const std::string & primer : *sets[s]){
                auto iter = pSet.find(primer);
                if(iter == pSet.end()){
                    continue;
                }
                int idx = std::distance(iter0,iter);
                char buf[100];
                sprintf(buf,"Primer_%04d",idx);
                if(!strandList[s].empty()){
                    strandList[s] += ",";
                }
                strandList[s] += buf;
                nPrimer[s]++;
            }
        }
        targetInfo << id << "\t" << nPrimer[0]*nPrimer[1] << "\t" << strandList[0] << "\t" << strandList[1] << "\n";
    }
    targetInfo.close();
    Log("Logged primer and target metadata",DEBUG);
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}


//Identify Manditory Primers (Those which target targets targeted by only one primer)
void id_mandatory_primers(PrimerSet & mSet,TargetSet & tSet,TargetByPrimer &p_tMap,PrimerByTarget &t_pMap){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    Log("Identifying manditory Primers ...",INFO);
    for(const auto & x : t_pMap){
        tSet.insert(x.first);
        if(x.second.size() == 1){
            mSet.insert(x.second[0]);
        }
    }
    //Remove any targets covered by the mandatory set from consideration
    for(const std::string & primer: mSet){
        for(const auto & target : p_tMap[primer]){
            tSet.erase(target.to_str());
        }
    }
    Log("Identifed %d manditory Primers ...",DEBUG,mSet.size());
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}

/*
void mapping_blacklist(std::string file, TargetByPrimer &p_tMap,PrimerByTarget &t_pMap){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    if(file.empty()){
        return;
    }
    Log("Excluding primers inappropriately mapping to the blacklist (%s) ...",INFO,file.c_str());
    PrimerSet pSet;
    for(const auto & x : p_tMap){
        pSet.insert(x.first);
    }
    std::string candidateFile = OutputPrefix + "_candidates.fna";
    std::string mapFile = OutputPrefix + "_blacklist.sam";
    std::ofstream candidateStream(candidateFile.c_str(),std::ios::out | std::ios::trunc);
    std::string bwaCommand = "bwa mem" +
        " -t " + std::to_string(Threads) + " -k " + std::to_string(MinLen / 2) +
        " -o " + mapFile + " " + file + " " + candidateFile;
    output_primer_set(pSet,candidateStream);
    candidateStream.close();
    if(system(bwaCommand.c_str())){
        Log("Error in mapping to blacklist",ERROR);
        exit(1);
    }
    seqan3::sequence_file_input fin{mapFile};
    using record_type = decltype(fin)::record_type;
    std::unordered_map<std::string,std::vector<BlacklistRegion> blMap{};
    for (auto & record : fin){
        if(record.flag() & 0x4){ //Skip unmapped reads
            continue;
            auto rID = record.reference_id();
            auto rPos = record.reference_position();
            BlacklistRegion 
            bool strand = (record.flag() & 0x10);
            auto seq = record.sequence();
        }
    }
    Log("Primers Blacklisted",DEBUG);
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}
*/

template<typename ... Args>
void Log(const char * format, LogLevel lvl, Args ... args){
#ifdef _TIMING
    auto clockStart = std::chrono::steady_clock::now();
#endif
    std::string header("[");
    char buffer[100];
    time_t rawTime;
    time(&rawTime);
    struct tm * timeinfo = localtime(&rawTime);
    strftime(buffer,100,"%x %R",timeinfo);
    header += buffer;
    header += "] ";
    switch(lvl){
        case ERROR:
            header += "ERROR";
            break;
        case WARNING:
            header += "WARNING";
            break;
        case INFO:
            header += "INFO";
            if(!bVerbose && !bDebug){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
                return;
            }
            break;
        case DEBUG:
            header += "DEBUG";
            if(!bDebug){
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
                return;
            }
            break;
        default:
            throw "Invalid LogLevel\n";
    }
    header += " - ";
    int size_s = std::snprintf( nullptr, 0, (header + format).c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, (header + format).c_str(), args ... );
    std::cerr << std::string(buf.get(), buf.get() + size - 1) << "\n";
#ifdef _TIMING
        auto clockEnd = std::chrono::steady_clock::now();
        ExecutionTime[__func__] += std::chrono::duration_cast<std::chrono::nanoseconds>(clockEnd - clockStart).count();
#endif
}




#ifdef _TEST
void run_test(){
    PrimerByTarget t_pMap{};
    TargetByPrimer p_tMap{};

    //Init t_pMap, and a set of targets
    std::vector<struct TargetInfo> alltargets{};
    for(int i = 0; i < 5; i++){
        struct TargetInfo target;
        target.id = std::to_string(i);
        target.strand = false;
        alltargets.push_back(target);
        target.strand=true;
        alltargets.push_back(target);
    }
    for(const struct TargetInfo & target : alltargets){
        //std::cerr << target.id << target.to_str() << "\n";
        t_pMap[target.to_str()] = std::vector<std::string>();
    }
    t_pMap[alltargets[0].to_str()].push_back(std::string("ACATAGAA"));
    t_pMap[alltargets[0].to_str()].push_back(std::string("E"));
    t_pMap[alltargets[0].to_str()].push_back(std::string("F"));

    t_pMap[alltargets[1].to_str()].push_back(std::string("ACATAGAA"));

    t_pMap[alltargets[2].to_str()].push_back(std::string("ACATAGAA"));
    t_pMap[alltargets[2].to_str()].push_back(std::string("E"));

    t_pMap[alltargets[3].to_str()].push_back(std::string("GACATGGC"));
    t_pMap[alltargets[3].to_str()].push_back(std::string("G"));
    t_pMap[alltargets[3].to_str()].push_back(std::string("H"));

    t_pMap[alltargets[4].to_str()].push_back(std::string("GACATGGC"));
    t_pMap[alltargets[4].to_str()].push_back(std::string("G"));
    t_pMap[alltargets[4].to_str()].push_back(std::string("H"));

    t_pMap[alltargets[5].to_str()].push_back(std::string("GACATGGC"));
    t_pMap[alltargets[5].to_str()].push_back(std::string("G"));
    t_pMap[alltargets[5].to_str()].push_back(std::string("H"));


    t_pMap[alltargets[6].to_str()].push_back(std::string("GACATGGC"));
    t_pMap[alltargets[6].to_str()].push_back(std::string("TTCTGGCGA"));

    t_pMap[alltargets[7].to_str()].push_back(std::string("D"));
    t_pMap[alltargets[7].to_str()].push_back(std::string("I"));
    t_pMap[alltargets[7].to_str()].push_back(std::string("J"));
    t_pMap[alltargets[7].to_str()].push_back(std::string("K"));

    t_pMap[alltargets[8].to_str()].push_back(std::string("D"));
    t_pMap[alltargets[8].to_str()].push_back(std::string("I"));
    t_pMap[alltargets[8].to_str()].push_back(std::string("J"));
    t_pMap[alltargets[8].to_str()].push_back(std::string("K"));

    t_pMap[alltargets[9].to_str()].push_back(std::string("D"));
    t_pMap[alltargets[9].to_str()].push_back(std::string("I"));
    t_pMap[alltargets[9].to_str()].push_back(std::string("J"));
    t_pMap[alltargets[9].to_str()].push_back(std::string("K"));

    //Init p_tMap
    for(std::string primer : {"ACATAGAA","GACATGGC","TTCTGGCGA","D","E","F","G","H","I","J","K"}){
        p_tMap[primer] = std::vector<struct TargetInfo>();
    }
    p_tMap[std::string("ACATAGAA")].push_back(alltargets[0]);
    p_tMap[std::string("ACATAGAA")].push_back(alltargets[1]);
    p_tMap[std::string("ACATAGAA")].push_back(alltargets[2]);

    p_tMap[std::string("GACATGGC")].push_back(alltargets[3]);
    p_tMap[std::string("GACATGGC")].push_back(alltargets[4]);
    p_tMap[std::string("GACATGGC")].push_back(alltargets[5]);
    p_tMap[std::string("GACATGGC")].push_back(alltargets[6]);

    p_tMap[std::string("TTCTGGCGA")].push_back(alltargets[6]);

    p_tMap[std::string("D")].push_back(alltargets[7]);
    p_tMap[std::string("D")].push_back(alltargets[8]);
    p_tMap[std::string("D")].push_back(alltargets[9]);

    p_tMap[std::string("E")].push_back(alltargets[0]);
    p_tMap[std::string("E")].push_back(alltargets[2]);

    p_tMap[std::string("F")].push_back(alltargets[0]);

    p_tMap[std::string("G")].push_back(alltargets[3]);
    p_tMap[std::string("G")].push_back(alltargets[4]);
    p_tMap[std::string("G")].push_back(alltargets[5]);

    p_tMap[std::string("H")].push_back(alltargets[3]);
    p_tMap[std::string("H")].push_back(alltargets[4]);
    p_tMap[std::string("H")].push_back(alltargets[5]);

    p_tMap[std::string("I")].push_back(alltargets[7]);
    p_tMap[std::string("I")].push_back(alltargets[8]);
    p_tMap[std::string("I")].push_back(alltargets[9]);

    p_tMap[std::string("J")].push_back(alltargets[7]);
    p_tMap[std::string("J")].push_back(alltargets[8]);
    p_tMap[std::string("J")].push_back(alltargets[9]);

    p_tMap[std::string("K")].push_back(alltargets[7]);
    p_tMap[std::string("K")].push_back(alltargets[8]);
    p_tMap[std::string("K")].push_back(alltargets[9]);

    init_subseq_gc(p_tMap);

    PrimerSet pSet{"ACATAGAA","GACATGGC","TTCTGGCGA"};    

    double curLoss = 0.0;
    update_loss_function(curLoss,PrimerSet(),PrimerSet(),pSet);
    std::cerr << curLoss << "\n";
    std::cerr << evaluate_loss_function(pSet) << "\n";
    pSet.erase("ACATAGAA");
    update_loss_function(curLoss,pSet, PrimerSet({"ACATAGAA"}), PrimerSet());
    std::cerr << curLoss << "\n";
    std::cerr << evaluate_loss_function(pSet) << "\n";
    pSet.clear();
    update_loss_function(curLoss,pSet,PrimerSet({"GACATGGC","TTCTGGCGA"}), PrimerSet({"GACATGGC"}));
    pSet.insert("GACATGGC");
    std::cerr << curLoss << "\n";
    std::cerr << evaluate_loss_function(pSet) << "\n";
}
#endif
