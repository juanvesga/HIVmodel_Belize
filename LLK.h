#pragma once
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>


#ifndef LLK_h
#define LLK_h

using namespace std;

const int nLLKs = 6;


// function declarations

double* LLK(double * inparams, int sizep, int time0, int time_end, double llks_out[]);
double unifrnd(double min, double max);
double llk_bino(double sample, double datapoint, double modelpoint);
double llk_poiss(double data, double model);
double llk_poiss2(double data, double model);
double llk_unif(double datapoint, double modelpoint, double lowDom, double upDom);
double llk_unif2(double datapoint, double modelpoint);
double llk_negbin(double datapoint, double modelpoint);
double llk_lognr(double datapoint, double modelpoint);
bool Bounds(std::vector<double> params, vector<vector<double>> &bounds);
double LogPrior(std::vector<double> parameters, string distribution, int dim, vector<vector<double>> &bounds, double a[] = 0, double b[] = 0);
void step_sizes(std::vector<vector<double>> &proposal_width, double step, vector<vector<double>> &bounds);
void random_params(std::vector<double> &param_array, double nparams, vector<vector<double>> &bounds);
double sqdiff(double simu, double data);
#endif /* LLK_h */


#pragma once
