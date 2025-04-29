#include "Snnt_constants.h"
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "Class/SNN.h"

static vector<double> PreSpike_Time;
static vector<int> PreSpike_Stream;
static vector<int> PreSpike_Signal; // 0 for background hit, 1 for signal hit, 2 for L1 neuron spike
static vector<int> PreSpike_Class;

using namespace std;


static int N_part; // Number of generated particles in an event
static double First_angle;
static int ievent;

// clear hits vector
void Reset_hits()
{
    hit_pos.clear();
    return;
}

// Start binning r-z plane ---------------
int GetBinR(double r_hit)
{
    if (r_hit < 0 || r_hit > max_R)
        return N_bin_r - 1;

    // 10 bins are associated to a tracking layer
    for (int i = 0; i < N_TrackingLayers; i++)
    {
        if (r_hit > Left_Layers[i] && r_hit < Right_Layers[i])
            return i;
        else if (r_hit < Left_Layers[i+1])
            return i;
    }
    return N_bin_r - 1;
}

int GetBinZ(double z)
{
    double tmp = (z + z_range / 2.);
    if (tmp < 0)
        tmp = 0;
    else if (tmp > z_range)
        tmp = z_range - epsilon;

    return (int)(tmp / z_bin_length);
}

int GetStreamID(int r, int z)
{
    return r + z * N_bin_r; // Visually streams are sorted by z and then by r
    // return r*N_bin_z+z;  // Visually streams are sorted by r and then by z
}
// End binning r-z plane -----------------

// Transforms hits into spikes streams by scanning the event
void Encode(double t_in)
{
    // sort by phi angle
    sort(hit_pos.begin(), hit_pos.end(), [](const Hit &h1, const Hit &h2)
         { return h1.phi < h2.phi; });

    for (auto &&row : hit_pos)
    {
        double time = t_in + row.phi / omega;
        int itl = GetStreamID(GetBinR(row.r), GetBinZ(row.z));

        PreSpike_Time.push_back(time);
        PreSpike_Stream.push_back(itl);
        PreSpike_Signal.push_back(row.id); //1,2 -> respectively Backgroung, Signal
        PreSpike_Class.push_back(row.pclass);
    }

    // rescan from [0, delta]
    for (auto &&row : hit_pos)
    {
        if (row.phi > delta)
            break;
        double time = t_in + (row.phi + M_PI * 2.) / omega;

        int itl = GetStreamID(GetBinR(row.r), GetBinZ(row.z));

        PreSpike_Time.push_back(time);
        PreSpike_Stream.push_back(itl);
        PreSpike_Signal.push_back(row.id); // 1,2 -> respectively Backgroung, Signal
        if(row.pclass >= 0) PreSpike_Class.push_back(row.pclass + N_classes); //ghost particle -> I indicate it with a shifted pclass
        else PreSpike_Class.push_back(row.pclass);
    }
}

tuple<TFile*, TTree*, TTree*, TTree*> readRootFile(const string &rootInput)
{
    // Open the ROOT file in read mode
    TFile *file = TFile::Open(rootInput.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        throw runtime_error("Error: Cannot open file " + rootInput);
    }

    // Access the required directories
    TDirectoryFile *dirIT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidIT"));
    TDirectoryFile *dirOT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidOT"));
    TDirectoryFile *dirEV = dynamic_cast<TDirectoryFile *>(file->Get("classification"));

    if (!dirIT)
    {
        file->Close();
        throw runtime_error("Error: Cannot access directory clusterValidIT in file " + rootInput);
    }

    if (!dirOT)
    {
        file->Close();
        throw runtime_error("Error: Cannot access directory clusterValidOT in file " + rootInput);
    }

    if (!dirEV)
    {
        file->Close();
        throw runtime_error("Error: Cannot access directory classification in file " + rootInput);
    }

    // Extract the trees from the directories
    TTree *IT = dynamic_cast<TTree *>(dirIT->Get("tree"));
    TTree *OT = dynamic_cast<TTree *>(dirOT->Get("tree"));
    TTree *ET = dynamic_cast<TTree *>(dirEV->Get("event_tree"));

    if (!IT)
    {
        file->Close();
        throw runtime_error("Error: Cannot access tree in clusterValidIT in file " + rootInput);
    }

    if (!OT)
    {
        file->Close();
        throw runtime_error("Error: Cannot access tree in clusterValidOT in file " + rootInput);
    }

    if (!ET)
    {
        file->Close();
        throw runtime_error("Error: Cannot access event_tree in classification in file " + rootInput);
    }

    // Do not close the file here; it will be managed by the caller
    return {file, IT, OT, ET};
}

// To read our preprocessed file
void ReadFromProcessed(TTree *IT, TTree *OT, TTree *ET, int id_event_value)
{
    Reset_hits();
    First_angle = max_angle;

    //retrieve the information about the eventClass
    int eventClass;
    ET->SetBranchAddress("eventClass", &eventClass);
    ET->GetEntry(id_event_value-1);
    pclass = eventClass;
    
    //TODO: generalize to more than 2 particles with an integer division
    if(eventClass < 0) {N_part=0; pclass=0;}  
    else if(eventClass >= 0 && eventClass < N_classes) N_part=1;
    else N_part = 2;
    
    // Reading the Inner Tracker

    float z;
    float r, phi;
    float id_event;
    float type;
    float cluster_pclass;
    
    IT->SetBranchAddress("cluster_z", &z);
    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("eventID", &id_event);
    IT->SetBranchAddress("cluster_type", &type);
    IT->SetBranchAddress("pclass", &cluster_pclass);

    if (ievent % NROOT == 0)
    {
        last_row_event_IT = 0;
        last_row_event_OT = 0;
    }

    // Loop over entries and find rows with the specified id_event value
    for (int i = last_row_event_IT; i < IT->GetEntries(); ++i)
    {
        IT->GetEntry(i);

        if (static_cast<int>(id_event) != id_event_value)
        {
            last_row_event_IT = i;
            break;
        }
        phi += M_PI;
        if (static_cast<int>(type) == 1)
        {
            type = SIG;
            phi += 2. * M_PI * ((int)(ievent / (NROOT))) * 1. / ((int)(N_events / NROOT) + 1);
            if (phi >= 2. * M_PI)
                phi -= 2. * M_PI;
            if (phi < First_angle)
                First_angle = phi;
        }
        else
        {
            type = BGR;
            if (phi >= 2. * M_PI)
                phi -= 2. * M_PI;
        }

        hit_pos.emplace_back(r, z, phi, static_cast<int>(type), cluster_pclass);
    }

    // OUT Tracker

    OT->SetBranchAddress("cluster_z", &z);
    OT->SetBranchAddress("cluster_R", &r);
    OT->SetBranchAddress("cluster_phi", &phi);
    OT->SetBranchAddress("eventID", &id_event);
    OT->SetBranchAddress("cluster_type", &type);
    OT->SetBranchAddress("pclass", &cluster_pclass);

    for (int i = last_row_event_OT; i < OT->GetEntries(); ++i)
    {
        OT->GetEntry(i);
        if (static_cast<int>(id_event) != id_event_value)
        {
            last_row_event_OT = i;
            break;
        }
        phi += M_PI;
        if (static_cast<int>(type) == 1)
        {
            phi += 2. * M_PI * ((int)(ievent / (NROOT))) * 1. / ((int)(N_events / NROOT) + 1);
            type = SIG;
            if (phi >= 2. * M_PI)
                phi -= 2. * M_PI;
            if (phi < First_angle)
                First_angle = phi;
        }
        else
        {
            type = BGR;
            if (phi >= 2. * M_PI)
                phi -= 2. * M_PI;
        }

        hit_pos.emplace_back(r, z, phi, static_cast<int>(type), cluster_pclass);
    }
}

vector<vector<int>> buildSpikeMatrix(const vector<int>& PreSpike_Stream, const vector<double>& PreSpike_Time) {
    // Get the size of the vectors
    int num_rows = N_TrackingLayers;
    int num_cols = (int) (max_angle/omega);
    cout << "num_cols: " << num_cols << endl;

    // Initialize the matrix with zeros
    vector<vector<int>> spike_matrix(num_rows, vector<int>(num_cols, 0));

    // Fill the matrix with ones at the specified positions
    for (int k = 0; k < PreSpike_Stream.size(); ++k) {
        int i = PreSpike_Stream[k]; // Row index
        int j = static_cast<int>(PreSpike_Time[k]); // Column index (cast to int)
        
        // Ensure the indices are within valid bounds
        if (i >= 0 && i < num_rows && j >= 0 && j < num_cols) {
            spike_matrix[i][j] = 1;
        } else {
            
            stringstream ss;
            ss << "Error: Indices are out of matrix bounds. "
               << "num_rows=" << num_rows << ", num_cols=" << num_cols 
               << ", i=" << i << ", j=" << j << ", k=" << k;
            throw out_of_range(ss.str());
        }
    }

    return spike_matrix;
}

void printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}


void writeMatrixToReadableFile(const vector<vector<int>>& matrix, const string& filename) {
    // Open the file in write mode
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }

    // Loop through each row and write it to the file
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            file << matrix[i][j];
            if (j < matrix[i].size() - 1) {
                file << ",";  // Add a comma between values
            }
        }
        file << "\n";  // Newline after each row
    }

    // Close the file after writing
    file.close();
}

void writeMatrixToFile(const vector<vector<int>>& matrix, const string& filename) {
    // Open the file in binary write mode
    ofstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }

    // Get the number of rows and columns from the matrix
    size_t num_rows = matrix.size();
    size_t num_cols = matrix[0].size();

    // Loop through each column and write the column data
    for (size_t col = 0; col < num_cols; ++col) {
        unsigned short val = 0;  // 16-bit value to store packed data for the column
        for (size_t row = 0; row < num_rows; ++row) {
            val = (val << 1) | matrix[row][col];  // Pack each bit into the 16-bit value (shift left and OR)
        }
        file.write(reinterpret_cast<char*>(&val), sizeof(val));  // Write the packed column data as 2 bytes
    }

    // Close the file after writing
    file.close();
}

void Produce_Stimuli(string rootInput="./DATA/muons_amuons_100k_100br.root", int _N_events=1, bool debugging=true){
    auto [file, IT, OT, ET] = readRootFile(rootInput);
    for(ievent=0; ievent < _N_events; ievent++){
        if (ievent % NROOT == 0)
        {
            last_row_event_IT = 0;
            last_row_event_OT = 0;
        }

        ReadFromProcessed(IT, OT, ET, ievent % NROOT + 1);

        double t_in = ievent * (max_angle + Empty_buffer) / omega; // Initial time -> every event adds 25 ns
        Encode(t_in);
        if (debugging)
        {
            cout << "Event " << ievent <<  "encoded" << endl;
            cout << "   Stream class time" << endl;
            for (int ih = 0; ih < PreSpike_Class.size(); ih++)
            {
                cout << "   " << PreSpike_Stream[ih] << " " << PreSpike_Class[ih] << " " << PreSpike_Time[ih] << endl; 
            }
        }
        try 
        {
            vector<vector<int>> spike_matrix = buildSpikeMatrix(PreSpike_Stream, PreSpike_Time);
            writeMatrixToReadableFile(spike_matrix, "matrix.csv");
            writeMatrixToFile(spike_matrix, "matrix.bin");
        } 
        catch (const out_of_range& e) 
        {
            cerr << "Error: " << e.what() << endl;
        }
    }
}
