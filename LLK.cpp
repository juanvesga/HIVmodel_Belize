#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <string>
#include "boost/random.hpp"
#include <random>
#include "boost/math/distributions/uniform.hpp"
#include "boost/math/distributions/binomial.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include "C:\Users\JFV09\Dropbox\TB\Code\TBmax_LHS\LLK.h"
#include "C:\Users\JFV09\Documents\Visual Studio 2015\Projects\TBmax_LHS\TBmax_LHS\make_modelsimple.h"
#include "C:\Users\JFV09\Documents\Visual Studio 2015\Projects\TBmax_LHS\TBmax_LHS\model3.h"


using namespace std;

double* LLK(double * inparams, int sizep, int time0, int time_end, double llks_out[]) {


double inparams[nParams] = { 0 };
			for (int j = 0; j < nParams; j++) {
				inparams[j] = LHStable[i][j];
			}


			//	// Make a model class
			Model * modelzero = new Model;
			//
			//	// Run model with initial set of parameters (setZero)
			Results * output_0 = modelzero->Run(inparams, nParams,time0, time_end , intvn);



			int checkisnan = 0, reject = 0;
			//******************************************************* Pass Model outputs into new arrays

			double msm_model[2] = { 0.0 };
			double gen_model[2] = { 0.0 };
			double preg_model[2] = { 0.0 };
			double art_model= 0.0;
			double plha_model =0;
			double fsw_model = 0;

			msm_model[0] = output_0->state[prevmsmall][yearsLimmsm[0] - timeZero];
			msm_model[1] = output_0->state[prevmsmall][yearsLimmsm[1] - timeZero];

			gen_model[0] = output_0->state[prev][yearsLimgen[0] - timeZero];
			gen_model[1] = output_0->state[prev][yearsLimgen[1] - timeZero];

			preg_model[0] = output_0->state[prevpreg49][yearsLimpreg[0] - timeZero];
			preg_model[1] = output_0->state[prevpreg49][yearsLimpreg[1] - timeZero];

			art_model = (output_0->state[onart][yearsLimART - timeZero]) + (output_0->state[suppr][yearsLimART - timeZero]) + (output_0->state[diseng][yearsLimART - timeZero]);
			plha_model = output_0->state[plha][yearsLimplha - timeZero];

			fsw_model = output_0->state[prevfsw][yearsLimfsw - timeZero];


			checkisnan = isnan(msm_model[0]) + isnan(fsw_model) + isnan(art_model)+ isnan(gen_model[0]) ;
						
			int checks = 0;
			if (msm_model[0] < 0.0 || art_model < 0.0 || gen_model[0]< 0.0 )
			{
				checks = checks + 1;
			}


			//if (mortisnan + prevmsmisnan + prevfswisnan + prevsisnan > 0) { cout << "NaN or INF" << endl; };
			if (checkisnan> 0) {
				reject = reject + 1; //.........
			}
			for (int j = 0; j < DurSim; j++) {
				checks += (output_0->state[inc][j] < 0);
				//if (checks > 0) { cout << "checkI" << endl; }
			}
	//******************************************************* Pass Model outputs into new arrays & likelihoods

			double prevgen_model[2] = { 0,0 };
			double prevmsm_model = 0;
			double prevfsw_model = 0;
			double prevpreg_model[5] = {0.0 };
			double hivcases_model[25]={0.0};
			double onART_model[4] = { 0.0 };
			double mort_model[9] = { 0.0 };

			for (int j = 0; j < onART_years.size(); j++) {
				int  year = 0;
				year = onART_years[j] - timeZero;
				onART_model[j] = (output_0->state[onart][year])  + (output_0->state[suppr][year])+ (output_0->state[diseng][year]);
			}
			int  year = msm_years[0] - timeZero;
				prevmsm_model = (output_0->state[prevmsmall][year]);

		        year = fsw_years[0] - timeZero;
				prevfsw_model = (output_0->state[prevfsw][year]);
			
			for (int j = 0; j < preg_years.size(); j++) {
				int  year = 0;
				year = preg_years[j] - timeZero;
				prevpreg_model[j] = (output_0->state[prevpreg49][year]);
			}

			for (int j = 0; j < hivcases_years.size(); j++) {
				int  year = 0;
				year = hivcases_years[j] - timeZero;
				hivcases_model[j] = output_0->state[nothiv_a0][year] + output_0->state[nothiv_a1][year] + output_0->state[nothiv_a2][year] + output_0->state[nothiv_a3][year] +
					output_0->state[notaids_a0][year] + output_0->state[notaids_a1][year] + output_0->state[notaids_a2][year] + output_0->state[notaids_a3][year];
			}
			for (int j = 0; j < mort_years.size(); j++) {
				int  year = 0;
				year = mort_years[j] - timeZero;
				mort_model[j] = (output_0->state[mortaAIDS][year]);
			}


		//     Clear-up heap....
	delete model;
	delete recent;

	// Reject if Bad Runs
	if ((out_isnan + out_isinf + out_isneg) > 0) {
		double fail[nLLKs] = { 0 };

		for (int g = 0; g<nLLKs; g++) {
			fail[g] = 2141454.0;

		}
		return fail;
	};


			/////////////////////////////////\\\\\\\\\\\\\ LIKELIHOODS ///////??????///////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

			double llkart = 0.0, llkmsm = 0.0, llkfsw = 0.0, llkpreg = 0.0;
			double llkcases = 0.0, llkmort=0.0;

			for (int j = 0; j < onART_years.size(); j++) {
					llkart += llk_poiss(onART_data[j], onART_model[j]); ////cout << 1 << endl;
			}
				llkmsm = llk_bino(130, prevmsm_data, prevmsm_model); 
				llkfsw = llk_poiss(prevfsw_data, prevfsw_model); 
			
			for (int j = 0; j < preg_years.size(); j++) {
				llkpreg += llk_poiss(prevpreg_data[j], prevpreg_model[j]); ////cout << 6 << endl;
			}
			for (int j = 0; j < hivcases_years.size(); j++) {
				llkcases += llk_poiss(hivcases_data[j], hivcases_model[j]);
			}
			for (int j = 0; j < mort_years.size(); j++) {
				llkmort += llk_poiss(mort_data[j], mort_model[j]);
			}






	llks_out[0] = llkmsm;
	llks_out[1] = llkfsw;
	llks_out[2] = llkpreg;
	llks_out[3] = llkart;
	llks_out[4] = llkcases;
	llks_out[5] = llkmort;
	

	return llks_out;
}

////////////////////Uniform Distributed Random numbers within  range//////////////////////////////////////////////
double unifrnd(double min, double max) {
	std::uniform_real_distribution<double> distribution(min, max);
	std::random_device rd;
	std::default_random_engine generator(rd());

	double response;

	response = distribution(generator);
	return response;
}

/////////////Create a vector with uniform starter paramns

void random_params(vector<double> &param_array, double nparams, vector<vector<double>> &bounds)
{
	for (int ii = 0; ii < nparams; ii++) {

		param_array[ii] = unifrnd(bounds[0][ii], bounds[1][ii]);


	}
}

/////////////////  ==========================================================//////////////////////////////////////////////
////////////LogLikelihood Binomial distribution
double llk_bino(double sample, double datapoint, double modelpoint) {
	double llkbino = 0.0;
	if (modelpoint < 0.0) { llkbino = -10000.0; std::cout << "negative model, sample " << sample << endl; }
	else {
		boost::math::binomial_distribution <> myBino(sample, datapoint);
		llkbino = log(pdf(myBino, modelpoint*sample));

		if (std::isinf(llkbino) == 1) { llkbino = -500.0; }
	}
	return llkbino;

}
// LogLike poisson
double llk_poiss(double data, double model) {
	double llkpois = 0.0;
	if (model < 0.0) { llkpois = -50.0; std::cout << "negative model poissART " << endl; }
	else {
		llkpois = log(exp(data * log(model) - lgamma(data + 1.0) - model));
		if (std::isinf(llkpois) == 1) { llkpois = -200.0; }
	}
	return llkpois;
}
// LogLike poisson
double llk_poiss2(double data, double model) {
	double llkpois = 0.0;
	if (model < 0.0) { llkpois = -1000.0;  std::cout << "negative model poissART " << endl; }
	else {
		llkpois = log(exp(data * log(model) - lgamma(data + 1.0) - model));
		if (std::isinf(llkpois) == 1) { llkpois = -1000.0; }
	}
	return llkpois;
}
// LogLike Uniform
double llk_unif(double datapoint, double modelpoint, double lowDom, double upDom) {
	double llkunif = 0.0;
	if (modelpoint < 0.0) { std::cout << "negative model" << endl; }
	if (modelpoint < 0.0) { llkunif = -15; }
	else {
		boost::math::uniform_distribution <> myUnif(lowDom, upDom);
		llkunif = log(pdf(myUnif, modelpoint));
		if (std::isinf(llkunif) == 1) { llkunif = -15.0; }
	}
	return llkunif;
}


double llk_unif2(double datapoint, double modelpoint) {
	double llkunif = 0.0;
	double lowDom = datapoint - (datapoint*0.01);
	double upDom = datapoint + (datapoint*0.01);
	if (modelpoint < 0.0) { std::cout << "negative model" << endl; }
	if (modelpoint < 0.0) { llkunif = -30; }
	else {
		boost::math::uniform_distribution <> myUnif(lowDom, upDom);
		llkunif = log(pdf(myUnif, modelpoint));
		if (std::isinf(llkunif) == 1) { llkunif = -30.0; }
	}
	return llkunif;
}
double llk_negbin(double datapoint, double modelpoint) {
	double llknegbin = 0.0;
	int point = floor(modelpoint);

	if (point < 0) { std::cout << "negative model, sample " << endl; }
	if (point < 0) { llknegbin = -800; std::cout << "negative model, sample " << endl; }
	else {
		//std::default_random_engine generator;
		boost::math::negative_binomial_distribution<> myNeg(datapoint, 0.5);
		llknegbin = log(pdf(myNeg, point));
		if (std::isinf(llknegbin) == 1) { llknegbin = -800.0; }
	}
	return llknegbin;

}
// Log Normal Distribution 
double llk_lognr(double datapoint, double modelpoint) {
	double llklognr = 0.0;
	if (modelpoint < 0.0) { llklognr = -10000.0; std::cout << "negative model " << endl; }
	else {
		boost::math::lognormal_distribution <> mylgnr(log(datapoint), 0.45);
		llklognr = log(pdf(mylgnr, modelpoint));

		if (std::isinf(llklognr) == 1) { llklognr = -500.0; }
	}

	return llklognr;
}


// file bounds
bool Bounds(vector<double> params, vector<vector<double>> &bounds) {
	int size = params.size();
	int violations = 0;
	bool result;
	int ids[] = { 4, 5 };
	for (int j = 0; j < size; j++) {
		if (params[j] < 0) {
			violations = violations + 1;
			//cout << "negative: " << j << endl;
		}
	}

	for (int jj = 0; jj < 2; jj++) {
		int j = ids[jj];
		if (!((params[j] <= bounds[1][j]) & (params[j] >= bounds[0][j]))) {
			violations = violations + 1;
			//	cout << "out of range: " << j<<" "<<params[j]<< endl;
		}
	}

	result = (violations > 0);

	return result;

}

// LogPriors

double LogPrior(vector<double> parameters, string distribution, int dim, vector<vector<double>> &bounds, double a[], double b[])

/*Function to calculate logPriors.
Inputs:
parameters : vector of parameter values
a: if distrbution normal a corresponds to means , if uniform is the left domain
b: if distrbution normal a corresponds to std , if uniform is the right domain
dimension: number of parameters
*/


{

	double logprior = 0;


	if (distribution == "Normal") {

		for (int j = 0; j < dim; j++) {

			double mean = a[j];
			double sd = b[j];
			boost::math::normal_distribution <> myN(mean, sd);
			double val = 0;
			val = log(pdf(myN, parameters[j]));
			if (isinf(val) == 1) {
				val = -1;
				logprior = logprior + val;
			}
			else {

				logprior = logprior + val;
			}
		}


	}



	if (distribution == "Uniform") {


		for (int j = 0; j < dim; j++) {
			double low = bounds[0][j];
			double up = bounds[1][j];
			boost::math::uniform_distribution <> myUnif(low, up);
			double val = 0;
			val = log(pdf(myUnif, parameters[j]));
			if (isinf(val) == 1) {
				val = -1;
				logprior = logprior + val;
			}
			else {

				logprior = logprior + val;
			}
		}



	}


	return logprior;



}


void step_sizes(vector<vector<double>> &proposal_width, double step, vector<vector<double>> &bounds)
{
	int dim = proposal_width[0].size();
	for (int j = 0; j < dim; j++) {
		//	for (int i = 0; i < dim; i++) {

		proposal_width[j][j] = (bounds[1][j] - bounds[0][j])*step;


		//	}
	}

}
double sqdiff(double simu, double data) {
	double ls = log((pow((data - simu), 2)));

	return ls;

}