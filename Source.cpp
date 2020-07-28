
/////////////////BLOCK WISE MCMC -MH Procedure with proposal adaptation using variance covariance matrix
// created by Juan F Vesga 13/04/2016

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <random>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include "boost/random.hpp"
#include "C:\Users\JFV09\documents\visual studio 2015\Projects\Belizemcmc\Belizemcmc\LLK.h"
#include "C:\Users\JFV09\documents\visual studio 2015\Projects\Belizemcmc\Belizemcmc\mvrnorm.h"
#include "file_path.h"
#pragma comment(linker, "/STACK:200000000")
#pragma comment(linker, "/HEAP:200000000")


using namespace std;


int main(int argc, char* argv[]) {
	///////////>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>test
	int nSamples = 0;
	int burn_in = 0;
	int adapt_l = 0;//Length of adaptaion phases
	int adapt_n = 0;// number of adaptation phases
	int	random_start = 0;
	int chainNo = 1;
	double step = 0;
	double init_scale = 1;
	//string region;
	if (argc > 1) {
		nSamples = atoi(argv[1]);
		burn_in = atoi(argv[2]);
		adapt_l = atoi(argv[3]);
		adapt_n = atoi(argv[4]);
		random_start = atoi(argv[5]);
		chainNo = atoi(argv[6]);
		step = atof(argv[7]);

	}
	else {
		nSamples = 1000;
		burn_in = 0;
		adapt_l = 0;//Length of adaptaion phases
		adapt_n = 0;// number of adaptation phases
		random_start = 0;
		chainNo = 0;
		step = 0.05;
	}

	const int nParams = 24;
	const int time0 = 1975;
	const int time_end = 2017;
	const int DurSim = 1 + (time_end - time0);
	double step_var = step;
	double scale = init_scale;
	int nchains = 3;
	int usecovmat = 0;

	//~~~~~~~~Parameter sampling boundaries
	// Bounds 
	std::vector<std::vector<double>>Bounds{
	{1E-6,0.002,0.05,0.01,0,0.00013,0.0006,2,0.1,2000,0.01,2,1,2,120,0,0,2,1,0.5,0.5,0.5,0.5,4.7},
	{1E-3,0.01,0.12,0.1,1, 0.00141, 0.00109,20,0.6,2017,0.5,20,10,20,300,4,4,20,10,1,1.5,1.5,1.5,18.8}	
    }; // beta_urb, beta_rur,  careseek_urb , careseek_rur,   RR_h 

	   //~~~~~~~~Pick  files and dataheader

	go_tofile();

				 //~~~~~~~~~~~~~~Read Files 
	double* sdZero = new double[nParams];
	double * startbest = new double[nParams];
	double * init_mean = new double[nParams];
	double** BestRuns = new double*[nParams];
	double** covMat = new double*[nParams];
	double* posteriors = new double[nSamples];
	double* acceptance = new double[nSamples];
	double* accRate = new double[nSamples];
	double** LLKs = new double*[nSamples];
	double** params = new double*[nSamples];
	double** propo = new double*[nSamples];
	double** sdvec = new double*[nSamples];
	

	for (int i = 0; i < nSamples; i++) {
		LLKs[i] = new double[nLLKs];
		params[i] = new double[nParams];
		propo[i] = new double[nParams];
		sdvec[i] = new double[nParams];
	
		posteriors[i] = 0;
		accRate[i] = 0;
		acceptance[i] = 0;

		for (int j = 0; j < nParams; j++) {
			params[i][j] = 0;
			propo[i][j] = 0;
			sdvec[i][j] = 0;

		}
		for (int j = 0; j < nLLKs; j++) {
			LLKs[i][j] = 0;
		}

	}
	for (int j = 0; j < nParams; j++) {
		covMat[j] = new double[nParams];
		for (int i = 0; i < nParams; i++) {

			covMat[j][i] = 0.0;
		}
	}
	for (int j = 0; j < nParams; j++) {
		sdZero[j] = 0.0;
		startbest[j] = 0.0;
		init_mean[j] = 0.0;
		BestRuns[j] = new double[nchains];

		for (int jj = 0; jj < nchains; jj++) {

			BestRuns[j][jj] = 0.0;

		}
	}
	//Vectors
	int vartrackLength = adapt_l; 	//VARIABLE LENGTHS OF ADAPTIVE PHASES

	vector<double> sample_mean(nParams);
	vector<vector<double>> sample_var;
	sample_var.resize(nParams, vector<double>(nParams, 0.0));
	vector<vector<double>> proposal_width;
	proposal_width.resize(nParams, vector<double>(nParams, 0.0));
	vector<double> sum_x(nParams);
	vector<double> sum_dx(nParams);
	vector<vector<double>> sum_dxy;
	sum_dxy.resize(nParams, vector<double>(nParams, 0.0));
	vector<double>  InParamValues(nParams, 0); //**initialze params for first run
	vector<double> param_current(nParams, 0);

	ifstream inputFile2;
	inputFile2.open(sigmafile[0], ios::out);
	if (!inputFile2)
		std::cout << "The file can't be opened";
	else {
		while (!inputFile2.eof())
		{
			for (int j = 0; j < nParams; j++) {
				inputFile2 >> sdZero[j];
				if (j == nParams - 1) { inputFile2.close(); } //close the file when it is complete read out.

			}
		}
	}

	ifstream bestFile;
	bestFile.open(bestrunfile[0], ios::out);
	if (!bestFile)
		std::cout << "The file can't be opened";

	else {

		while (!bestFile.eof())

		{
			for (int k = 0; k < nParams; k++) {
				for (int j = 0; j < nchains; j++) {
					bestFile >> BestRuns[k][j];
				}
				if (k == nParams - 1) { bestFile.close(); } //close the file when it is complete read out.
			}

		}
	}


	ifstream meanfile;
	meanfile.open(mean_file[0], ios::out);
	if (!meanfile)
		std::cout << "The file can't be opened";

	else {

		while (!meanfile.eof())

		{
			for (int k = 0; k < nParams; k++) {

				meanfile >> init_mean[k];
				if (k == nParams - 1) { meanfile.close(); } //close the file when it is complete read out.
			}

		}
	}

	ifstream covFile;
	covFile.open(covmat_file[0], ios::out);
	if (!covFile)
		std::cout << "The file can't be opened";

	else {

		while (!covFile.eof())

		{
			for (int k = 0; k < nParams; k++) {
				for (int j = 0; j < nParams; j++) {
					covFile >> covMat[k][j];
				}
				if (k == nParams - 1) { covFile.close(); } //close the file when it is complete read out.
			}

		}
	}


	// Use covariance matrix from previous runs
	if (usecovmat == 1) {
		for (int j = 0; j < nParams; j++) {
			for (int i = 0; i < nParams; i++) {
				//covMat[j][i] = covMat[j][i] * pscale;
				proposal_width[i][j] = covMat[i][j];
			}
		}
	}
	else {
		step_sizes(proposal_width, step_var, Bounds);

	}


	if (random_start == 1) {
		vector<double> param_array(nParams, 0.0);

		random_params(param_array, nParams, Bounds);

		for (int j = 0; j < nParams; j++) {
			startbest[j] = param_array[j];
		}
	}
	else {

		for (int j = 0; j < nParams; j++) {
			startbest[j] = BestRuns[j][chainNo];
		}

	}

	// ********************* Prepare for MCMC

	for (int h = 0; h < nParams; h++) {

		InParamValues[h] = startbest[h];
	}
	for (int i = 0; i < nParams; i++)  //initialise values used for aptaptive phases
	{
		sample_mean[i] = 0.0;
		sum_dx[i] = 0.0;
		sum_x[i] = 0;
		for (int j = 0; j < nParams; j++)
		{
			sample_var[i][j] = 0.0;
			sum_dxy[i][j] = 0.0;
		}
	}

	for (int i = 0; i < nParams; i++)  //initialise covariance matrix with the diagonal of starting SD
	{
		sdZero[i] = (proposal_width[i][i]); //sdZero[i] * sdfac;
	}
	for (int k = 0; k < nParams; k++) {
		params[0][k] = InParamValues[k];
	}
	int n = 0;
	int n_adapt = 1;
	double acc = 0;
	double acc_lag = 0;

	// **********************************************************  Initial Run

	//*****          get LogPriors
	double globalLogPrior = LogPrior(InParamValues, "Uniform", nParams, Bounds);

	//*******Call likelihood fucntion (data | model)
	double  out_LLK[nLLKs] = { 0.0 };
	double* loglikelihoods = LLK(startbest, nParams, time0, time_end, out_LLK);

	double globalLLK = 0;
	for (int i = 0; i < nLLKs; i++) {
		globalLLK += loglikelihoods[i];
	}
	double currentPosterior = globalLogPrior + globalLLK;

	posteriors[0] = currentPosterior;
	LLKs[0][0] = loglikelihoods[0];// msm;
	LLKs[0][1] = loglikelihoods[1];//fswr;
	LLKs[0][2] = loglikelihoods[2];// preg;
	LLKs[0][3] = loglikelihoods[3];//art;
	LLKs[0][4] = loglikelihoods[4];//cases;
	LLKs[0][5] = loglikelihoods[5];//Mort;
	
	param_current = InParamValues;
	//***************************************** Start Iterations *****************************************
	for (int l = 0; l < nSamples - 1; l++)
	{
		accRate[l] = acc / double(l);	if (isnan(accRate[l]) | isinf(accRate[l])) { accRate[l] = 0; }
		///////////////////MH ALgorithm////////////////////////////////////////////
		vector<double> proposals(nParams, 0);

		// Re-set values of params and posterior for next iteration
		posteriors[l + 1] = currentPosterior;
		for (int k = 0; k < nParams; k++) {
			params[l + 1][k] = params[l][k];

		}
		for (int k = 0; k < nLLKs; k++) {
			LLKs[l + 1][k] = LLKs[l][k];
		}

		//// Propose new set of parameters
		//if (l <= burn_in + vartrackLength | adapt_n==0) {// Before adaptation do simple paramater drawing
		//	for (int j = 0; j < nParams; j++) {
		//		proposals[j] = params[l][j] + normrnd(0, 1)*proposal_width[j][j]*scale;
		//		}
		//}
		//else { // adapt and after
		proposals = multinormal_sample(nParams, proposal_width, param_current, scale);

		//		}
		//cout << param_current[2] / pscale << endl;
		//cout << proposals[2]/pscale << endl;
		for (int j = 0; j < nParams; j++) {
			propo[l][j] = proposals[j];
		}

		// check proposals inside Bounds
		if (Boundsfun(proposals, Bounds)) {
			//proposals = param_current;
			//cout << l << " " << "\nOut of Bounds" << ":\n";
			continue;
		}

		// get proposals vector into array to fit LLK definition (Lazy)
		double* proposals_arr = &proposals[0];


		//get LogPriors
		globalLogPrior = LogPrior(proposals, "Uniform", nParams, Bounds);


		//Call likelihood fucntion (data | model)
		double  out_LLK[nLLKs] = { 0.0 };
		double* loglikelihoods_new = LLK(proposals_arr, nParams, time0, time_end, out_LLK);

		// Check if return erros& continue
		if (loglikelihoods_new[0] == 2141454) {
			//cout << "Failed" << endl;
			continue;
		}

		double globalLLK = 0;
		for (int i = 0; i < nLLKs; i++) {
			globalLLK += loglikelihoods_new[i];
		}
		double proposedPosterior = globalLogPrior + globalLLK;
		//check errors and don't let it pass
		if (globalLLK >= 0 || globalLLK < -1000000) proposedPosterior = posteriors[l];


		//Acceptance Procedure
		double log_acceptance = 0.0;
		double r = 0.0;
		log_acceptance = proposedPosterior - currentPosterior;
		if (log_acceptance == 0) { log_acceptance = -100; } // Prevent parameters out of bound to increase acc rate
		r = log(unifrnd(0.0, 1.0));
		if (r <= log_acceptance) {		// if accepted:
			for (int j = 0; j < nParams; j++) {
				params[l + 1][j] = proposals[j];
			}
			param_current = proposals;
			currentPosterior = proposedPosterior;
			posteriors[l + 1] = proposedPosterior;
			LLKs[l + 1][0] = loglikelihoods_new[0];// msm;
			LLKs[l + 1][1] = loglikelihoods_new[1];// fsw;
			LLKs[l + 1][2] = loglikelihoods_new[2];// preg;
			LLKs[l + 1][3] = loglikelihoods_new[3];// ART;
			LLKs[l + 1][4] = loglikelihoods_new[4];// cases;
			LLKs[l + 1][5] = loglikelihoods_new[5];// Mort;

			acc++;
		}//******************* END of MH ALGORITHM


		 //////////////////////////////// Proposal adaptation using MVRNORM

		if (l > burn_in && l < burn_in + adapt_n*vartrackLength)
		{
			n++;
			var_track(&param_current, &sample_mean, &sample_var, &sum_x, &sum_dxy, nParams, n);
		}
		if (l > burn_in && l < burn_in + adapt_n*vartrackLength)
		{
			n++;
			var_track(&param_current, &sample_mean, &sample_var, &sum_x, &sum_dxy, nParams, n);
		}
		if (l == burn_in + n_adapt*vartrackLength && n_adapt <= adapt_n)
		{

			double acceptanceRate = accRate[l];// (acc - acc_lag) / double(n);
											   /////Corey Chivers adapt of MVRN
			coVar(&sample_var, &sum_dxy, nParams, n);
			if (is_pos_def(nParams, sample_var))
			{
				proposal_width = sample_var;
				cout << "\nProposal covariance matrix update " << n_adapt << ":\n";

			}
			else {
				cout << "\nNot Possitive Definate " << n_adapt << ":\n";
			}
			////*********************************

			/////// Alternatively, Re-scale k parameter for tuning//////////////////////

			/*if (!(acceptanceRate >= 0.15 && acceptanceRate <= 0.5)) {

			tune_scale(scale, acceptanceRate);

			}*/
			/////////////////////////////////////////////////////

			n_adapt++;
			//re-initialise values used for adaptive algorithm
			n = 0;
			acc_lag = acc;
			for (int j = 0; j < nParams; j++)
			{
				sample_mean[j] = param_current[j];
				sum_dx[j] = 0;
				sum_x[j] = 0;
				for (int k = 0; k < nParams; k++)
					sum_dxy[j][k] = 0;
			}

		}



	}//*********************End of iterations **************************/



	 //fill last row
	for (int j = 0; j < nParams; j++) {
		acceptance[nSamples - 1] = acceptance[(nSamples - 2)];
		accRate[nSamples - 1] = accRate[(nSamples - 2)];
	}


	// Write results to text
	fstream file;
	if (argc > 1) { file.open(argv[8], ios::out); }
	//    else { file.open("/Users/juan/Dropbox/TB/Code/mcmc/mcmcDebug.txt", ios::out); }//mac
	else { file.open("mcmcDebug.txt", ios::out); }//pc

	for (int t = 0; t < nSamples; t++) {
		file << params[t][0] << " "
			<< params[t][1] << " "
			<< params[t][2] << " "
			<< params[t][3] << " "
			<< params[t][4] << " "
			<< params[t][5] << " "
			<< params[t][6] << " "
			<< params[t][7] << " "
			<< params[t][8] << " "
			<< params[t][9] << " "
			<< params[t][10] << " "
			<< params[t][11] << " "
			<< params[t][12] << " "
			<< params[t][13] << " "
			<< params[t][14] << " "
			<< params[t][15] << " "
			<< params[t][16] << " "
			<< params[t][17] << " "
			<< params[t][18] << " "
			<< params[t][19] << " "
			<< params[t][20] << " "
			<< params[t][21] << " "
			<< params[t][22] << " "
			<< params[t][23] << " "
			<< sdvec[t][0] << " "
			<< sdvec[t][1] << " "
			<< sdvec[t][2] << " "
			<< sdvec[t][3] << " "
			<< sdvec[t][4] << " "
			<< sdvec[t][5] << " "
            << sdvec[t][6] << " "
			<< sdvec[t][7] << " "
			<< sdvec[t][8] << " "
			<< sdvec[t][9] << " "
			<< sdvec[t][10] << " "
			<< sdvec[t][11] << " "
             << sdvec[t][12] << " "
			<< sdvec[t][13] << " "
			<< sdvec[t][14] << " "
			<< sdvec[t][15] << " "
			<< sdvec[t][16] << " "
			<< sdvec[t][17] << " "
             << sdvec[t][18] << " "
			<< sdvec[t][19] << " "
			<< sdvec[t][20] << " "
			<< sdvec[t][21] << " "
			<< sdvec[t][22] << " "
			<< sdvec[t][23] << " "
			<< propo[t][0] << " "
			<< propo[t][1] << " "
			<< propo[t][2] << " "
			<< propo[t][3] << " "
			<< propo[t][4] << " "
			<< propo[t][5] << " "
			<< propo[t][6] << " "
			<< propo[t][7] << " "
			<< propo[t][8] << " "
			<< propo[t][9] << " "
			<< propo[t][10] << " "
			<< propo[t][11] << " "
			<< propo[t][12] << " "
			<< propo[t][13] << " "
			<< propo[t][14] << " "
			<< propo[t][15] << " "
			<< propo[t][16] << " "
			<< propo[t][17] << " "
			<< propo[t][18] << " "
			<< propo[t][19] << " "
			<< propo[t][20] << " "
			<< propo[t][21] << " "
			<< propo[t][22] << " "
			<< propo[t][23] << " "
			<< accRate[t] << " "
			<< posteriors[t] << " "
			<< LLKs[t][0] << " "
			<< LLKs[t][1] << " "
			<< LLKs[t][2] << " "
			<< LLKs[t][3] << " "
			<< LLKs[t][4] << " "
			<< LLKs[t][5] << " " << endl;
	}
	file.close();

	return 0;
}
