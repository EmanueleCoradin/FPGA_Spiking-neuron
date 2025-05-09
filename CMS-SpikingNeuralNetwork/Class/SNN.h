#ifndef SNN_H
#define SNN_H

#include <iostream>
#include <algorithm>

#include <fstream>
#include <cstdlib>
#include <string>

#include <vector>
#include <stdexcept>
#include <functional>
#include "TRandom3.h"
#include "nlohmann/json.hpp"

using namespace std;

class SNN
{
private:
    double bisectionMethod(double a, double b, int in, double epsilon, std::function<float(int, double, bool)> func);

public:
    const int EPSP  = 1;
    const int IPSP  = 2;
    const int SPIKE = 0;
    const int NOCLASS  = -2;
    const int BKGCLASS = -1;
    const int SIGCLASS = 0;
    float alpha;              // 0.25; // factor tuning inhibition strength
    
    float CFI0;
    float CFI1;
    float CF01;

    float L1inhibitfactor;      // multiplier for L1 inhibits

    float K;
    float K1;                   // constants to tune post-synaptic shape
    float K2;                   // see above

    double IPSP_dt_dilation;     // shape factor of exponential IPSP signal

    double MaxDelay;         // Determines shape of IE signal[s]

    double tau_m;          // membrane time constant[s]
    double tau_s;          // synapse time constant[s]
    double tau_r;          // refractory time[s]
    double tau_plus;       // [s]
    double tau_minus;      // [s]
    double taud_plus;      // [s]
    double taud_minus;     // [s]
    double taud_plus_2;    // [s]
    double taud_minus_2;   // [s]

    double a_plus;          
    double a_minus;       
    double d_plus;         
    double d_minus;       

    float sparsity;
    bool split_layer0;      //decide if you want to split layer 0 in two.

    int N_InputStreams;
    int N_streams;

    double fire_granularity;    //it defines how much close we will look for neuron's activation.
    float fire_precision;      // [V] it defines the precision of the neuron firetime detection.

    double Delta_delay;
    double Mean_delay;

    //Variables that depend on the upper ones
    int N_neurons;
    float Threshold[2];     // Neuron threshold in arbitrary units; in paper it is 550V but they have 1000 channels, 100x ours
    double tmax;
    double MaxDeltaT;    // time window wherein pre-synaptic, inhibition, and post-synaptic kernels affect neuron potential


    /* data */
    //------- variables ---------
    float **Weight;         // Weight of synapse-neuron strength
    float **Weight_initial; // Initial configuration of the weights
    bool **check_LTD;       // checks to generate LTD after neuron discharge
    bool **Void_weight;     // These may be used to model disconnections
    double **Delay;          // Delay in incoming signals
    double **Delay_initial;
    bool **EnableIPSP;      // N_Neurons * N_Neurons matrix to turn off IPSP inside layers, "splitting" them in more sectors.

    vector<double> *History_time;       // Time of signal events per each 1neuron
    vector<int> *History_type;         // Type of signal
    vector<int> *History_ID;           // ID of generating signal stream or neuron
    vector<pair<int, int>> *History_ev_class;
    
    vector<double> *Fire_time;          // Times of firing of each neuron
    int *Neuron_layer;
    float *sumweight; // summed weights of streams for each neurons for the purpose of normalization
    double *sumdelays; // summed delays of streams for each neurons for the purpose of normalization
    
    int N_neuronsL[2];           // Number of neurons in layers 0 and 1
    TRandom3 *myRNG;
    double largenumber;
    double epsilon;
    double time_scaleFactor;

    SNN();
    SNN( int _NL0, int _NL1,
         float _alpha,
         float _CFI0, float _CFI1, float _CF01,
         float _L1inhibitfactor,
         float _K, float _K1, float _K2,
         double _IPSP_dt_dilation,
         double _MaxDelay,

         double _tau_m, double _tau_s, double _tau_r, double _tau_plus, double _tau_minus,
         double _a_plus, double _a_minus,

         double _taud_plus, double _taud_minus,
         double _taud_plus_2, double _taud_minus_2,
         double _d_plus, double _d_minus,

         int _N_InputStreams,
         float _Threshold0, float _Threshold1, float _sparsity, bool _split_layer0);
         
    ~SNN();


    //------- Functions ---------
    void Reset_Parameters(float _alpha,
         float _L1inhibitfactor,
         float _K, float _K1, float _K2,
         double _IPSP_dt_dilation,
         double _MaxDelay,

         double _tau_m, double _tau_s, double _tau_r, double _tau_plus, double _tau_minus,
         double _a_plus, double _a_minus,

         double _taud_plus, double _taud_minus,
         double _taud_plus_2, double _taud_minus_2,
         double _d_plus, double _d_minus,

         float _Threshold0, float _Threshold1);
    void Init_neurons(int ievent);
    void Init_weights_uniform();
    void Init_weights();
    void Set_weights();
    void Reset_weights();

    void Init_delays_man();
    void Init_delays_gauss();
    void Init_delays_uniform();

    void insert_spike(int id_neuron, double spike_time, int type, int id, int spike_class, int ievent);
    
    void Init_connection_map();
    float EPS_potential(double delta_t);
    float Spike_potential(double delta_t, int ilayer);
    float Inhibitory_potential(double delta_t, int ilayer);
    double Neuron_firetime(int in, double t);
    vector<pair <int, int>> Inspect_History(int in, double fire_time, double window);
    void Activate_Neuron(int in, double t);
    float Neuron_Potential(int in, double t, bool delete_history);
    void LTP_weights(int in, double fire_time, bool nearest_spike_approx, SNN &old);  
    void LTP_delays(int in, double fire_time, bool nearest_spike_approx, SNN &old);  
    void LTD_weights(int in, double fire_time, bool nearest_spike_approx, SNN &old);
    void LTD_delays(int in, double fire_time, bool nearest_spike_approx, SNN &old);
    void Renorm_weights(int in, SNN &old);
    void Renorm_delays(int in, SNN &old);
    void PrintWeights();
    void PrintDelays();
    void PrintSNN();
    void copy_from(const SNN& other);
    void dumpToJson(const string& filename);
    void loadFromJson(const string& filename);

    void Exclude_neuron(int in);
    void Exclude_neurons(const vector<int>& neurons_to_exclude);
    void CopyNeuronFrom(int from, int to);
    void MirrorDelays(int in);
};
#endif