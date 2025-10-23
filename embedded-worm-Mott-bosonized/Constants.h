#pragma once
#include "SmoWoSampler.h"  //Choose from "SmoWoSampler.h", "WoSampler.h"


//Physics
std::vector<int> Ls = { 256 };
std::vector<int> Bs = { 256 };
std::vector<double> Ks = { 0.35 };
std::vector<double> mus = { 0 };
double g = 1;

//Sampling
double lambda_w = 1;
double lambda_r_prefactor = 0.1;
std::string path = "C:\\Users\\oscar\\algo_Mott\\results\\SmoWo_tests";    //for instance C:\\Users\\...
std::vector<std::string> OBSERVABLES({ "N_space", "N_time", "field" });  //Select any from: "N_space", "N_time", "kappa", "rho_s", "C_2kf", "C_theta", "Ct_phi", "Cx_phi", "algotime", "field"
int PERCENT_SAVE = 1;
double TOTAL_NUMBER_SAMPLE = 1000;
double OUTPUT_DISTANCE_SWEEPS = 20;
int data_block_avg = 20;