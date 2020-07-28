//
//  modelTBsear_cpp
//  TBsear
//
//  Created by juan fernando vesga on 01/03/2017.
//  Copyright © 2017 juan fernando vesga. All rights reserved.
// This model uses an age structure divided by age ranges 


#include "Agefile.hpp"
#include "BirthRate.h"
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include "mortalityratesF.hpp"
#include "mortalityratesM.hpp"
#include "model.hpp"
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;
using namespace boost::numeric::odeint;


////----------------------- Type definitions-------------------------------------------------------------------------------------------
typedef boost::array< double, S_all > state_type;
typedef vector<double> input;
//--------------------------------------------------------------------------------------------------------------------------------------

////----------------------- Auxiliary Functions----------------------------------------------------------------------------------------
// interpolation
double interp(double t, double x0, double x1, double y0, double y1) {
	double output = y0 + (y1 - y0) * (t - x0) / (x1 - x0);
	return output;
}

// Swith index for x(RxSxA)
int switch_id(int s, int r, int a) {
	int index = s * R * A + r * A + a;
	return index;
}
// Differential after integration (for cuulative stages
vector<double> diff(vector<double> data)
{
	vector<double> vect_diff;

	for (unsigned int i = 0; i < data.size(); i++)
	{
		if (i == 0) {
			vect_diff.push_back(data[0]);
		}
		else {
			vect_diff.push_back(data[i] - data[i - 1]);
		}
	}
	return vect_diff;
}

// Results array creator

Results::Results(int Length) :
	arraySize(Length)
{
	CreateArray();

}
Results::~Results()

{}
//define custom function
void Results::CreateArray()
{
	//Create array
	int width = arraySize;
	int height = S_all;
	state.resize(height);
	for (int i = 0; i < height; i++) state[i].resize(width, 0.0);
	// return the array
	return;
}
//--------------------------------------------------------------------------------------------------------------------------------------
//------------------------Parameters Function 
void regionParsgo(int reg_int) {

}

////----------------------- Ordinary Differential Equations-----------------------------------------------------------------------------
class odes
{
	input prm;
	//
public:
	odes(input G) : prm(G) {}
	void operator () (const state_type x, state_type &dxdt, const double t) {
		//-------------------------------------- Parameters
		double FSW_frac = prm[1]; //0.0005;//
		double CFSW_frac = prm[2]; //0.0005;//
		double MSM_frac = prm[3]; //0.0005;//
		double epsi = prm[4];
		double betaFtoM = prm[5];
		double betaMtoF = prm[6];
		double betaMSM = prm[6]*prm[7];
		double detectionrpryr = prm[8];
		double test_peak = prm[9];
		double ART_coef = prm[10];
		double c_msm = prm[11];
		double c_bi =  prm[12];
		double c_cfsw =  prm[13];
		double c_fsw =  prm[14];
		double c_m = prm[15];
		double c_f = prm[16];
		double actsMSM = prm[17];
		double condom_x = prm[18];
		double propBiacts =  prm[19];
		double condom_ymsm = prm[20];
		double condom_ysw = prm[21];
		double condom_yhet = prm[22];
		double Acute_coef = prm[23];
		//---------------------------------------------------------- :: Regional Specific parameters :: 
		double migrationrate;
		double *HIVratio = NULL;
		double *Dis500p = NULL;
		double *Dis500 = NULL;
		double *Dis350 = NULL;
		double *Dis200 = NULL;
		double *Back500p = NULL;
		double *Back500 = NULL;
		double *Back350 = NULL;
		double *Back200 = NULL;
		double *ART350p = NULL;
		double *ART350 = NULL;
		double *ART200 = NULL;
		double *Supp500p = NULL;
		double *Supp500 = NULL;
		double *Supp350 = NULL;
		double *Supp200 = NULL;
		double *proptoART = NULL;

		////>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.......................... SOUTH
		migrationrate = 0.85;// from 0 to inf (reflects positive or negative net migration. Factors birth rate)
		ART350p = new double[43]{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.4,0.43,0.46,0.48, 0.50, 0.50 };
		ART350 = new double[43]{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.01,0.02,0.05,0.05,0.08,0.1,0.12,0.15,0.25,0.29,0.33,0.36,0.4,0.43,0.46,0.48, 0.50, 0.50 };
		ART200 = new double[43]{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.01,0.02,0.05,0.05,0.08,0.1,0.12,0.15,0.25,0.29,0.33,0.36,0.4,0.43,0.46,0.48, 0.50, 0.50 };
		proptoART = new double[43]{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0.005255335,0.01051067,0.015766005,0.02102134,0.025227978,0.026779556,0.034465349,0.06116855,0.096124139,0.129415359,0.145070034,0.152981721,0.145471666,0.14604282,0.152030115,0.170637141,0.192427764,0.213133445,0.233760017,0.254504127 ,0.254504127 ,0.254504127 };

		// -------------------------------------Intervention Switches
		bool C = 1; // turns condoms On when == 1
		bool ARV = 1; //turns ART On when == 1
		bool notransm = 0; // 1= set transmission effect of ART to zero (normally set to 0)
		bool T = 1; //cuts path to tested
		bool nep = 1;
		bool nohiv = 1;// Turns HIV infection to zero : 0=no infection 
		bool pgrowth = 1; // turns pop growth to 0
		bool equi = 0; // turns model demography to equilibrium (1==equilibrium)
		//-------------------------------------- Detection Rate
		// Method 1 : cumulative logistic distribution 
		double test_scale = 0;
		double test = 0;
		double hivT = 0;
		double hivT_scale = 0;
		double mu_test = (test_start + test_peak) / 2; // Mean
		double e = exp(-abs(1.81379936423421785*(t - mu_test) / 2));
		if (t < test_start) {
			test_scale = 0;
		}
		else {
			if (t >= mu_test) {
				test_scale = 1. / (1. + e);
			}
			else {
				test_scale = e / (1 + e);
			}
		}
		if (test_scale > 0) {
			test = detectionrpryr*test_scale*T;
		}

		mu_test = (hivtest_start + test_peak) / 2; // Mean
		e = exp(-abs(1.81379936423421785*(t - mu_test) / 2));
		if (t < hivtest_start) {
			hivT = 0;
		}
		else {
			if (t >= mu_test) {
				hivT = 1. / (1. + e);
			}
			else {
				hivT = e / (1 + e);
			}
		}
		if (test_scale > 0) {
			test = detectionrpryr*test_scale*T;
			//cout <<t<<" "<< test <<" "<< endl;
		}
		if (hivT > 0) {
			hivT_scale = detectionrpryr*hivT*T;
			//cout <<t<<" "<< hivT_scale <<" "<< endl;
		}

		// --------------------------------------Natural History, mortality, progression

		double proghorz[] = { 1.0 / durAc, (dur500p[0] + dur500p[1]) / 2, (dur500[0] + dur500[1]) / 2, (dur350[0] + dur350[1]) / 2, (dur200[0] + dur200[1]) / 2 };
		double progvert[4];
		progvert[0] = 1 / (1 / proghorz[1] + (1 / proghorz[2]) + (1 / proghorz[3]) + (1 / proghorz[4]));
		progvert[1] = 1 / (1 / proghorz[2] + 1 / proghorz[3] + 1 / proghorz[4]);
		progvert[2] = 1 / (1 / proghorz[3] + 1 / proghorz[4]);
		progvert[3] = 1 / (1 / proghorz[4]);
		double progvertA[4];
		progvertA[0] = progvert[0] / surv_ARText;
		progvertA[1] = progvert[1] / surv_ARText;
		progvertA[2] = progvert[2] / surv_ARText;
		progvertA[3] = progvert[3] / surv_ARText;

		double mortSupp[4];
		mortSupp[0] = progvertA[0];
		mortSupp[1] = progvertA[1];
		mortSupp[2] = progvertA[2];
		mortSupp[3] = progvertA[3];


		//---------------------------------------------- ART, suppressed and disenga	ged Rates
		int i = int(t);// integer of t to used as index in pre-def vectors.
		int ii = int(t); if (t > 2017) { ii = 2016; }
		int iii = int(t) + 1; if (t > 2017) { iii = 2016; }
		double x0 = i;
		double x1 = i + 1;

		double treatment[4] = { 0.0 };// ART Rates 
		double p_art = 0.0; // proportion straight from detection to ART

//		double y350pi = 0, y350pii = 0, y350i = 0, y350ii = 0,y200i = 0, y200ii = 0;

		if (ii >= (1996) && ii < (2002)) {
			int f = 2002 - timeZero;
			treatment[0] = ART350p[f] * ARV;
			treatment[1] = ART350p[f] * ARV;
			treatment[2] = ART350[f] * ARV;
			treatment[3] = ART200[f] * ARV;
		}
		else
		{
			treatment[0] = interp(t, x0, x1, ART350p[ii - timeZero], ART350p[iii - timeZero]) *ARV;
			treatment[1] = interp(t, x0, x1, ART350p[ii - timeZero], ART350p[iii - timeZero])*ARV;
			treatment[2] = interp(t, x0, x1, ART350[ii - timeZero], ART350[iii - timeZero])*ARV;
			treatment[3] = interp(t, x0, x1, ART200[ii - timeZero], ART200[iii - timeZero])*ARV;
		}
		p_art = interp(t, x0, x1, proptoART[ii - timeZero], proptoART[iii - timeZero])*ARV;

		double back_art[4] = { 0.0 };// Back to ART Rates 
		back_art[0] = treatment[3] * backRR;
		back_art[1] = treatment[3] * backRR;
		back_art[2] = treatment[3] * backRR;
		back_art[3] = treatment[3] * backRR;

		double mu[R][A] = { {0.0} }, muF[4] = { 0.0 }, muM[4] = { 0.0 }, muFi[4] = { 0.0 }, muMi[4] = { 0.0 };
		for (int jj = 0; jj < 85; jj++) {
			if (jj < 10) {
				muF[0] += mortRatesF[jj][ii - timeZero];
				muM[0] += mortRatesM[jj][ii - timeZero];
				muFi[0] += mortRatesF[jj][iii - timeZero];
				muMi[0] += mortRatesM[jj][iii - timeZero];
			}
			if ((jj >= 10) & (jj < 35)) {
				muF[1] += mortRatesF[jj][ii - timeZero];
				muM[1] += mortRatesM[jj][ii - timeZero];
				muFi[1] += mortRatesF[jj][iii - timeZero];
				muMi[1] += mortRatesM[jj][iii - timeZero];
			}
			if ((jj >= 35) & (jj < 50)) {
				muF[2] += mortRatesF[jj][ii - timeZero];
				muM[2] += mortRatesM[jj][ii - timeZero];
				muFi[2] += mortRatesF[jj][iii - timeZero];
				muMi[2] += mortRatesM[jj][iii - timeZero];
			}
			if (jj >= 50) {
				muF[3] += mortRatesF[jj][ii - timeZero];
				muFi[3] += mortRatesF[jj][iii - timeZero];
				muM[3] += mortRatesM[jj][ii - timeZero];
				muMi[3] += mortRatesM[jj][iii - timeZero];
			}

		}

		//	cout<< mortRatesM[51][0] <<" "<<muMi[3]<<" " << muM[3]<<" "<< muFi[3]<<" " << muF[3]<<" "<< endl;

		muF[0] = muF[0] / 10;
		muM[0] = muM[0] / 10;
		muF[1] = muF[1] / 25;
		muM[1] = muM[1] / 25;
		muF[2] = muF[2] / 15;
		muM[2] = muM[2] / 15;
		muF[3] = muF[3] / 35;
		muM[3] = muM[3] / 35;

		muFi[0] = muFi[0] / 10;
		muMi[0] = muMi[0] / 10;
		muFi[1] = muFi[1] / 25;
		muMi[1] = muMi[1] / 25;
		muFi[2] = muFi[2] / 15;
		muMi[2] = muMi[2] / 15;
		muFi[3] = muFi[3] / 35;
		muMi[3] = muMi[3] / 35;

		for (int a = 0; a < A; a++) {
			mu[0][a] = interp(t, x0, x1, muF[a], muFi[a]);
			mu[1][a] = interp(t, x0, x1, muF[a], muFi[a]);
			mu[2][a] = interp(t, x0, x1, muF[a], muFi[a]);
			mu[3][a] = interp(t, x0, x1, muM[a], muMi[a]);
			mu[4][a] = interp(t, x0, x1, muM[a], muMi[a]);
			mu[5][a] = interp(t, x0, x1, muM[a], muMi[a]);
			mu[6][a] = interp(t, x0, x1, muM[a], muMi[a]);
			mu[7][a] = interp(t, x0, x1, muM[a], muMi[a]);

		}

		double birthrt = 0.0;
		birthrt = interp(t, x0, x1, birthrate[ii - timeZero], birthrate[iii - timeZero]);

		//--------------------------------------Condoms, Safe Syringe 

		double condmsm = fmin(cond_msm*condom_ymsm, 1);
		double condbi = fmin(cond_msm*condom_ymsm, 1);
		double condhet = fmin(cond_het*condom_yhet, 1);
		double condfsw = fmin(cond_fsw*condom_ysw, 1);

		double msm_scale = 0;
		msm_scale = fmin(fmax(t - cond_start, condmsm*0.25) / condom_x, condmsm);
		double condommsm = 0;
		condommsm = msm_scale;

		double bi_scale = 0;
		bi_scale = fmin(fmax(t - cond_start, condbi*0.25) / condom_x, condbi);
		double condombi = 0;
		condombi = bi_scale;


		double sw_scale = 0;
		sw_scale = fmin(fmax(t - cond_start, condfsw*0.25) / condom_x, condfsw);
		double condomsw = 0;
		condomsw = sw_scale;

		double het_scale = 0;
		het_scale = fmin(fmax(t - cond_start, condhet*0.25) / condom_x, condhet);
		double condomhet = 0;
		condomhet = het_scale;


		//-------------------------------------------- Partnership formation and mixing

		// Who Mixes with Who matrix _ heterosxual
		double wMw[R][R] = { { 0.0 } }, wMwM[R][R] = { { 0.0 } }, assort[R][R] = { { 0.0 } }, assortM[R][R] = { { 0.0 } };
		assort[lo_w][lo_m] = 1; //  	
		assort[fsw][cfsw] = 1;//
		assort[lo_m][lo_w] = 1;//
		assort[cfsw][fsw] = 1;//  
		assort[bi][lo_w] = 1;//
		assortM[bi][bi] = 1;
		assortM[msm][msm] = 1;

		wMw[lo_w][lo_m] = 1;
		wMw[lo_w][cfsw] = 1;
		wMw[lo_w][bi] = 1;
		wMw[fsw][lo_m] = 1;
		wMw[fsw][cfsw] = 1;
		wMw[fsw][bi] = 1;
		wMw[lo_m][lo_w] = 1;
		wMw[lo_m][fsw] = 1;
		wMw[cfsw][lo_w] = 1;
		wMw[cfsw][fsw] = 1;
		wMw[bi][lo_w] = 1;
		wMw[bi][fsw] = 1;
		//Homosexual
		wMwM[bi][bi] = 1;
		wMwM[bi][msm] = 1;
		wMwM[msm][bi] = 1;
		wMwM[msm][msm] = 1;

		double Prt_Ch[R] = { 0.0 };
		double Prt_ChM[R] = { 0.0 };
		Prt_Ch[lo_w] = c_f;   //women non FSW
		Prt_Ch[fsw] = c_fsw;   //FSW
		Prt_Ch[lo_m] = c_m;   //Men
		Prt_Ch[cfsw] = c_cfsw;   //
		Prt_Ch[bi] = c_bi*(1 - propBiacts);   //Bi
		Prt_ChM[bi] = c_bi*propBiacts;   //msm Low
		Prt_ChM[msm] = c_msm;   //msm High

		double pChange[R][A] = { { 0.0 } }, pChangeM[R][A] = { { 0.0 } };
		double Prt_age = 0.0;
		for (int a = 0; a < A; a++) {
			if (a == 0) { Prt_age = ratio15_25; }
			if (a == 1) { Prt_age = ratio25_50; }
			if (a == 2) { Prt_age = ratio50_65; }
			if (a == 3) { Prt_age = ratio65_99; }
			pChange[lo_w][a] = Prt_Ch[lo_w] * Prt_age;
			pChange[fsw][a] = Prt_Ch[fsw] * Prt_age;
			pChange[lo_m][a] = Prt_Ch[lo_m] * Prt_age;
			pChange[cfsw][a] = Prt_Ch[cfsw] * Prt_age;
			pChange[bi][a] = Prt_Ch[bi] * Prt_age;
			pChangeM[bi][a] = Prt_ChM[bi] * Prt_age;
			pChangeM[msm][a] = Prt_ChM[msm] * Prt_age;
		}

		//-------------------------Sex Acts-----------------------------//
		vector<int> male_ids = { lo_m,cfsw,bi };
		vector<int> female_ids = { lo_w,fsw };
		vector<int> nonsw_ids = { lo_m,bi };
		vector<int> nonsw_idsF = { lo_w };
		vector<int> msm_ids = { bi,msm };
		double acts[R][R] = { { 0.0 } };
				
		acts[lo_w][lo_m] = wMw[lo_w][lo_m] * actsHet;
		acts[lo_w][cfsw] = wMw[lo_w][cfsw] * actsHet;
		acts[lo_w][bi] =   wMw[lo_w][bi] * actsHet;
		acts[fsw][lo_m] =  wMw[fsw][lo_m] * actsHet;
		acts[fsw][cfsw] =  wMw[fsw][cfsw] * actsSW;
		acts[fsw][bi] =    wMw[fsw][bi] * actsSW;
		acts[lo_m][lo_w] = wMw[lo_m][lo_w] * actsHet;
		acts[lo_m][fsw] =  wMw[lo_m][fsw] * actsHet;
		acts[cfsw][lo_w] = wMw[cfsw][lo_w] * actsHet;
		acts[cfsw][fsw] =  wMw[cfsw][fsw] * actsSW;
		acts[bi][lo_w] =   wMw[bi][lo_w] * actsHet;
		acts[bi][fsw] =    wMw[bi][fsw] * actsSW;
		//Homosexual
		acts[bi][bi] = wMwM[bi][bi] * actsMSM;
		acts[bi][msm] = wMwM[bi][msm] * actsMSM;
		acts[msm][bi] = wMwM[msm][bi] * actsMSM;
		acts[msm][msm] = wMwM[msm][msm] * actsMSM;





		//---------------------Beta coefficients 
		double beta_co[S];
		for (int h = 0; h < S; h++) {
			beta_co[h] = nohiv*1.0;
		}
		beta_co[Ac] = nohiv*(Acute_coef); //acute
		beta_co[U_200] = nohiv*(Acute_coef*lateStageHR);// Late stage infectiousness
		beta_co[D_200] = nohiv*(Acute_coef*lateStageHR);// Late stage infectiousness
		beta_co[Di_200] = nohiv*(Acute_coef*lateStageHR);// Late stage infectiousness
		beta_co[A_500p] = nohiv*(ART_coef*0.5 + (notransm* (1 - ART_coef*0.5))); //Not suppresed VL
		beta_co[A_500] = nohiv*(ART_coef*0.5 + (notransm* (1 - ART_coef*0.5))); //Not suppresed VL
		beta_co[A_350] = nohiv*(ART_coef*0.5 + (notransm* (1 - ART_coef*0.5))); //Not suppresed VL
		beta_co[A_200] = nohiv*(ART_coef*0.5 + (notransm* (1 - ART_coef*0.5))); //Not suppresed VL
		beta_co[S_500p] = nohiv*(ART_coef + (notransm* (1 - ART_coef))); //Not suppresed VL
		beta_co[S_500] = nohiv*(ART_coef + (notransm* (1 - ART_coef))); //Not suppresed VL
		beta_co[S_350] = nohiv*(ART_coef + (notransm* (1 - ART_coef))); //Not suppresed VL
		beta_co[S_200] = nohiv*(ART_coef + (notransm* (1 - ART_coef))); //Not suppresed VL

		//----------------------------- Mixing matrix
		double xsum[R][A] = { {0.0} };
		double totPop = 0.0;
		for (int s = 0; s < S; s++) {
			for (int r = 0; r < R; r++) {
				for (int a = 0; a < A; a++) {
					int id = switch_id(s, r, a);
					xsum[r][a] += x[id];
					totPop += x[id];

				}
			}
		}



		double Tot_partnersRate[R][A] = { {0.0} };
		double Tot_partnersRateM[R][A] = { {0.0} };

		for (int a = 0; a < A; a++) {
			Tot_partnersRate[lo_w][a] = xsum[lo_w][a] * pChange[lo_w][a];//
			Tot_partnersRate[fsw][a] = xsum[fsw][a] * pChange[fsw][a];//
			Tot_partnersRate[lo_m][a] = xsum[lo_m][a] * pChange[lo_m][a];
			Tot_partnersRate[cfsw][a] = xsum[cfsw][a] * pChange[cfsw][a];
			Tot_partnersRate[bi][a] = xsum[bi][a] * pChange[bi][a];
			Tot_partnersRateM[bi][a] = xsum[bi][a] * pChangeM[bi][a];
			Tot_partnersRateM[msm][a] = xsum[msm][a] * pChangeM[msm][a];
		}

		// This matrix reflects the Sum of all possible partenrs *partnerchange rate
		// for every risk group. (is not the same denominator for every risk group..thats why)
		//Calculate Total No Partners
		double Sumpartners[R] = { 0.0 }, SumpartnersM[R] = { 0.0 };
		for (int r = 0; r < R; r++) {
			for (int a = 0; a < A; a++) {
				Sumpartners[lo_w] += Tot_partnersRate[r][a] * wMw[lo_w][r];
				Sumpartners[fsw] += Tot_partnersRate[r][a] * wMw[fsw][r];
				Sumpartners[lo_m] += Tot_partnersRate[r][a] * wMw[lo_m][r];
				Sumpartners[cfsw] += Tot_partnersRate[r][a] * wMw[cfsw][r];
				Sumpartners[bi] += Tot_partnersRate[r][a] * wMw[bi][r];
				SumpartnersM[bi] += Tot_partnersRateM[r][a] * wMwM[bi][r];
				SumpartnersM[msm] += Tot_partnersRateM[r][a] * wMwM[msm][r];
			}
		}

		double MixMatrixM[R][A][R][A] = { {{{0.0}}} };
		double MixMatrixF[R][A][R][A] = { {{{0.0}}} };
		double MixMatrixMSM[R][A][R][A] = { {{{0.0}}} };
		for (int kk = 0; kk < male_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < female_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int k = male_ids[kk];
						int j = female_ids[jj];
						MixMatrixM[k][a][j][h] = agepdf_short[a][h] * epsi*assort[k][j] +
							(1 - epsi)*(wMw[k][j] * Tot_partnersRate[j][h]) / Sumpartners[k];
					}
				}
			}
		}
		for (int kk = 0; kk < female_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < male_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int j = male_ids[jj];
						int k = female_ids[kk];
						MixMatrixF[k][a][j][h] = agepdfF_short[a][h] * epsi*assort[k][j] + (1 - epsi)* (wMw[k][j] * Tot_partnersRate[j][h]) / Sumpartners[k];
					}
				}
			}
		}
		for (int kk = 0; kk < msm_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < msm_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int k = msm_ids[kk];
						int j = msm_ids[jj];
						MixMatrixMSM[k][a][j][h] = agepdf_short[a][h] * epsi*assortM[k][j] +
							(1 - epsi)*(wMwM[k][j] * Tot_partnersRateM[j][h]) / SumpartnersM[k];
					}
				}
			}
		}
		//---------------------------------------------------- BALANCE MIXING MATRICES
		double delta[R][A][R][A] = { {{{0.0}}} };
		double pChnew[R][A][R][A] = { {{{0.0}}} };
		double deltaM[R][A][R][A] = { {{{0.0}}} };
		double pChnewMSM[R][A][R][A] = { {{{0.0}}} };
		for (int kk = 0; kk < male_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < female_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int k = male_ids[kk];
						int j = female_ids[jj];
						delta[k][a][j][h] = wMw[k][j] * ((MixMatrixF[j][h][k][a] * Tot_partnersRate[j][h]) / (MixMatrixM[k][a][j][h] * Tot_partnersRate[k][a]));

						if (isnan(delta[k][a][j][h]) == 1) {
							delta[k][a][j][h] = 0.0;
						}
						if (isinf(delta[k][a][j][h]) == 1) {
							delta[k][a][j][h] = 0.0;
						}
					}
				}
			}
		}
		for (int kk = 0; kk < male_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < female_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int k = male_ids[kk];
						int j = female_ids[jj];
						pChnew[k][a][j][h] = pChange[k][a] * pow(delta[k][a][j][h], theta[j]);
						pChnew[j][h][k][a] = pChange[j][h] * pow(delta[k][a][j][h], (-(1 - theta[j])));
						if (isnan(pChnew[k][a][j][h]) == 1) {
							pChnew[k][a][j][h] = 0.0;
						}
						if (isnan(pChnew[j][h][k][a]) == 1) {
							pChnew[j][h][k][a] = 0.0;
						}
						if (isinf(pChnew[k][a][j][h]) == 1) {
							pChnew[k][a][j][h] = 0.0;
						}
						if (isinf(pChnew[j][h][k][a]) == 1) {
							pChnew[j][h][k][a] = 0.0;
						}
					}
				}
			}
		}
		for (int kk = 0; kk < msm_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < msm_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int k = msm_ids[kk];
						int j = msm_ids[jj];
						deltaM[k][a][j][h] = wMwM[k][j] *
							((MixMatrixMSM[j][h][k][a] * pChangeM[j][h] * xsum[j][h]) /
								(MixMatrixMSM[k][a][j][h] * pChangeM[k][a] * xsum[k][a]));
					}
				}
			}
		}

		for (int kk = 0; kk < msm_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < msm_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int k = msm_ids[kk];
						int j = msm_ids[jj];
						pChnewMSM[k][a][j][h] = pChangeM[k][a] * pow(deltaM[k][a][j][h], theta[j]);
						pChnewMSM[j][h][k][a] = pChangeM[j][h] * pow(deltaM[k][a][j][h], (-(1 - theta[j])));
						if (isnan(pChnewMSM[k][a][j][h]) == 1) {
							pChnewMSM[k][a][j][h] = 0.0;
						}
						if (isnan(pChnewMSM[j][h][k][a]) == 1) {
							pChnewMSM[j][h][k][a] = 0.0;
						}
						if (isinf(pChnewMSM[k][a][j][h]) == 1) {
							pChnewMSM[k][a][j][h] = 0.0;
						}
						if (isinf(pChnewMSM[j][h][k][a]) == 1) {
							pChnewMSM[j][h][k][a] = 0.0;
						}

					}
				}
			}
		}
		// ------Condom use array RxR
		double condomuse[R][R] = { {0.0} };
		for (int jj = 0; jj < male_ids.size(); jj++) {
			int j = male_ids[jj];
			condomuse[lo_w][j] = condomhet * wMw[lo_w][j] * C;
			condomuse[fsw][j] = condomsw  * wMw[fsw][j] * C;

			condomuse[j][lo_w] = condomhet * wMw[j][lo_w] * C;
			condomuse[j][fsw] = condomsw * wMw[j][fsw] * C;
		}
		condomuse[msm][msm] = condommsm * wMwM[msm][msm] * C;
		condomuse[msm][bi] = condommsm * wMwM[msm][bi] * C;
		condomuse[bi][msm] = condommsm * wMwM[bi][msm] * C;
		condomuse[bi][bi] = condombi * wMwM[bi][bi] * C;

		//---------------------------------Force of Infection 
		//Transmission probability per stage of infection Female to male
		double popProp[R][S][A] = { {{0.0}} };
		double TransmissionRateInPartnerShipM[R][A][R][A] = { {{{0.0}}} };
		double TransmissionRateInPartnerShipF[R][A][R][A] = { {{{0.0}}} };
		double TransmissionRateInPartnerShipMSM[R][A][R][A] = { {{{0.0}}} };
		double FOI_HetermaleMat[R][A][R][A] = { {{{0.0}}} };
		double FOI_HeterfemaleMat[R][A][R][A] = { {{{0.0}}} };
		double FOI_MSMMat[R][A][R][A] = { {{{0.0}}} };
		double FOI_HetM_summed[R][A] = { {0.0} };
		double FOI_HetF_summed[R][A] = { {0.0} };
		double FOI_MSM_summed[R][A] = { {0.0} };
		double FOI_grand[R][A] = { {0.0} };
		for (int r = 0; r < R; r++) {
			for (int s = 1; s < S; s++) {
				for (int a = 0; a < A; a++) {
					int id = switch_id(s, r, a);
					popProp[r][s][a] = x[id] / xsum[r][a];
					if (isnan(popProp[r][s][a])) {
						popProp[r][s][a] = 0;
					}

				}
			}
		}

		// FOI Males
		for (int kk = 0; kk < male_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < female_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						for (int s = 1; s < S; s++) {
							int k = male_ids[kk];
							int j = female_ids[jj];
							TransmissionRateInPartnerShipM[k][a][j][h] += 1 - pow(1 - (popProp[j][s][h] * betaFtoM*beta_co[s] * (1 - condomuse[k][j] * Condom_eff)), acts[k][j]);
						}
					}
				}
			}
		}

		// FOI Females
		for (int kk = 0; kk < female_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < male_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						for (int s = 1; s < S; s++) {
							int j = male_ids[jj];
							int k = female_ids[kk];
							TransmissionRateInPartnerShipF[k][a][j][h] += 1 - pow(1 - (popProp[j][s][h] * betaMtoF*beta_co[s] * (1 - condomuse[k][j] * Condom_eff)), acts[k][j]);
						}

					}
				}
			}
		}


		// FOI MSM
		for (int kk = 0; kk < msm_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < msm_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						for (int s = 1; s < S; s++) {
							int k = msm_ids[kk];
							int j = msm_ids[jj];
							TransmissionRateInPartnerShipMSM[k][a][j][h] += 1 - pow(1 - (popProp[j][s][h] * betaMSM * beta_co[s] * (1 - condomuse[k][j] * Condom_eff)), acts[k][j]);

						}

					}
				}
			}
		}

		for (int kk = 0; kk < male_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < female_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int k = male_ids[kk];
						int j = female_ids[jj];
						FOI_HetermaleMat[k][a][j][h] = MixMatrixM[k][a][j][h] * pChnew[k][a][j][h] * TransmissionRateInPartnerShipM[k][a][j][h];
						FOI_HeterfemaleMat[j][h][k][a] = MixMatrixF[j][h][k][a] * pChnew[j][h][k][a] * TransmissionRateInPartnerShipF[j][h][k][a];
					}
				}
			}
		}
		for (int kk = 0; kk < msm_ids.size(); kk++) {
			for (int a = 0; a < A; a++) {
				for (int jj = 0; jj < msm_ids.size(); jj++) {
					for (int h = 0; h < A; h++) {
						int k = msm_ids[kk];
						int j = msm_ids[jj];
						FOI_MSMMat[k][a][j][h] = MixMatrixMSM[k][a][j][h] * pChnewMSM[k][a][j][h] * TransmissionRateInPartnerShipMSM[k][a][j][h];
					}
				}
			}
		}

		for (int k = 0; k < R; k++) {
			for (int a = 0; a < A; a++) {
				for (int j = 0; j < R; j++) {
					for (int h = 0; h < A; h++) {
						FOI_HetM_summed[k][a] += FOI_HetermaleMat[k][a][j][h];
						FOI_HetF_summed[k][a] += FOI_HeterfemaleMat[k][a][j][h];
						FOI_MSM_summed[k][a] += FOI_MSMMat[k][a][j][h];
					}
				}
			}
		}
		for (int k = 0; k < R; k++) {
			for (int a = 0; a < A; a++) {
				FOI_grand[k][a] = FOI_HetF_summed[k][a] + FOI_HetM_summed[k][a] + FOI_MSM_summed[k][a];
			}
		}

		//----------------------------- Maintain a fraction of indivuduals perptually undetected
		double undF = 0.0;
		double num = 0.0;
		double denom = 0.0;
		//Set arrays of K,S,Ac to zero
		for (int k = 0; k < R; k++) {
			for (int s = 1; s < S; s++) {
				for (int a = 0; a < A; a++) {
					int id = switch_id(s, k, a);
					if (s < D_500p) {
						num += x[id];
					}
					denom += x[id];

				}
			}
		}
		undF = num / denom;

		if (undF < Undetected_Fraction) {
			test = test * 0.9;
		}


		//---------------- Recruit rates, mortality , Aging 

		double next_Age[R*S*A] = { 0.0 };
		double Background_mortality[R] = { 0.0 };
		double deathNathist[R] = { 0.0 };
		double deathTested[R] = { 0.0 };
		double deathonART[R] = { 0.0 };
		double deathNathist200[R] = { 0.0 };
		double deathonART200[R] = { 0.0 };
		double deathTested200[R] = { 0.0 };
		double deathSupp[R] = { 0.0 };
		double deathSupp200[R] = { 0.0 };
		double deathDis[R] = { 0.0 };
		double deathDis200[R] = { 0.0 };
		double HIVmortality[R] = { 0.0 };
		double AIDSdeathonly[R] = { 0.0 };
		double TotalDeath[R] = { 0.0 };
		double Recruit[R] = { 0.0 };
		double N_group[R] = { 0.0 };
		for (int k = 0; k < R; k++) {
			for (int a = 0; a < A; a++) {
				Background_mortality[k] += xsum[k][a] * mu[k][a];
				int u1 = switch_id(U_500p, k, a); int u2 = switch_id(U_500, k, a); int u3 = switch_id(U_350, k, a); int u4 = switch_id(U_200, k, a);
				deathNathist[k] += x[u1] * mortHIV[0] + x[u2] * mortHIV[1] + x[u3] * mortHIV[2];
				deathNathist200[k] += x[u4] * mortHIV[3];
				int d1 = switch_id(D_500p, k, a); int d2 = switch_id(D_500, k, a); int d3 = switch_id(D_350, k, a); int d4 = switch_id(D_200, k, a);
				deathTested[k] += x[d1] * mortHIV[0] + x[d2] * mortHIV[1] + x[d3] * mortHIV[2];
				deathTested200[k] += x[d4] * mortHIV[3];
				int a1 = switch_id(A_500p, k, a); int a2 = switch_id(A_500, k, a); int a3 = switch_id(A_350, k, a); int a4 = switch_id(A_200, k, a);
				deathonART[k] += x[a1] * mortSupp[0] + x[a2] * mortSupp[1] + x[a3] * mortSupp[2];
				deathonART200[k] += x[a4] * mortSupp[3];
				int s1 = switch_id(S_500p, k, a); int s2 = switch_id(S_500, k, a); int s3 = switch_id(S_350, k, a); int s4 = switch_id(S_200, k, a);
				deathSupp[k] += (x[s1] * mortSupp[0]) + (x[s2] * mortSupp[1]) + (x[s3] * mortSupp[2]);
				deathSupp200[k] += x[s4] * mortSupp[3];
				int di1 = switch_id(Di_500p, k, a); int di2 = switch_id(Di_500, k, a); int di3 = switch_id(Di_350, k, a); int di4 = switch_id(Di_200, k, a);
				deathDis[k] += x[di1] * mortHIV[0] + x[di2] * mortHIV[1] + x[di3] * mortHIV[2];
				deathDis200[k] += x[di4] * mortHIV[3];
			}
		}
		double Alldead = 0.0;
		for (int k = 0; k < R; k++) {
			HIVmortality[k] = deathNathist[k] + deathTested[k] + deathonART[k] + deathSupp[k] + deathDis[k] + deathNathist200[k] + deathTested200[k] + deathonART200[k] + deathSupp200[k] + deathDis200[k];
			TotalDeath[k] = HIVmortality[k] + Background_mortality[k];
			Alldead += TotalDeath[k];
			AIDSdeathonly[k] = (deathTested200[k] + deathonART200[k] + deathSupp200[k] + deathDis200[k]) + deathNathist200[k] * undeR;
		}
		double propDead[R] = { 0.0 };
		for (int k = 0; k < R; k++) {
			propDead[k] = TotalDeath[k] / Alldead;

		}
		double births = 0.0, migrate = 0.0, input = 0.0;
		births = totPop * birthrt* migrationrate;

		//cout << totPop << endl;
		input = births;// +migrate - Alldead;
		double Q[R] = { 0.0 };
		Q[nr_w] = NoRW_frac * f_prop;//proportion of no risk females
		Q[lo_w] = (1 - (FSW_frac + NoRW_frac))*f_prop;// proportion Non-FSW Women
		Q[fsw] = FSW_frac * f_prop;// proportion of female sex Workers

		Q[nr_m] = NoRM_frac * m_prop;// proportion of no risk males
		Q[lo_m] = (1 - (NoRM_frac + CFSW_frac + MSM_frac))*m_prop; // proportion of heterosexual men
		Q[cfsw] = CFSW_frac * m_prop;// proportion of clients of female sex Workers
		Q[bi] = Bi_frac*MSM_frac * m_prop; // proportion of low risk MSM
		Q[msm] = (1 - Bi_frac)*MSM_frac * m_prop; // proportion of high risk MSM

		Recruit[0] = (TotalDeath[0] * equi) + Q[0] * input*(1 - equi);
		Recruit[1] = (TotalDeath[1] * equi) + Q[1] * input*(1 - equi);
		Recruit[2] = (TotalDeath[2] * equi) + Q[2] * input*(1 - equi);
		Recruit[3] = (TotalDeath[3] * equi) + Q[3] * input*(1 - equi);
		Recruit[4] = (TotalDeath[4] * equi) + Q[4] * input*(1 - equi);
		Recruit[5] = (TotalDeath[5] * equi) + Q[5] * input*(1 - equi);
		Recruit[6] = (TotalDeath[6] * equi) + Q[6] * input*(1 - equi);
		Recruit[7] = (TotalDeath[7] * equi) + Q[7] * input*(1 - equi);
		// ------------Population distribution
		double totpopage[A] = { 0 };
		for (int k = 0; k < R; k++) {
			for (int a = 0; a < A; a++) {
				totpopage[a] += xsum[k][a];
			}
		}
		//----------------------------Aging every year 	
		double modageRef = 1.0;
		double year_agecheck = fmod(t, modageRef);
		if (year_agecheck == 0.0) {

			std::memset(next_Age, 0, sizeof(next_Age[0])*R*S*A);


			for (int k = 0; k < R; k++) {
				for (int s = 0; s < S; s++) {
					int id1 = switch_id(s, k, 0);
					int id2 = switch_id(s, k, 1);
					int id3 = switch_id(s, k, 2);
					int id4 = switch_id(s, k, 3);
					next_Age[id1] = x[id1] - x[id1] / 10;
					next_Age[id2] = x[id2] + x[id1] / 10 - x[id2] / 25;
					next_Age[id3] = x[id3] + x[id2] / 25 - x[id3] / 15;
					next_Age[id4] = x[id4] + x[id3] / 35;

				}
			}

			for (int k = 0; k < R; k++) {
				for (int s = 0; s < S; s++) {
					for (int a = 0; a < A; a++) {
						int id = switch_id(s, k, a);
						dxdt[id] = next_Age[id];
					}
				}
			}
		}

		//--------------------------------------TURN OVER
		double turn_hr = 0.0;
		double swleaving = 0.0;
		double fhr = 0.0;
		double turnover_in[R][S][A] = { {{0.0 }} };
		double turnover_ou[R][S][A] = { {{0.0}} };
		double turnover[R][S][A] = { {{0.0}} };
		for (int s = 0; s < S; s++) {
			for (int a = 0; a < A; a++) {
				int id_fsw = switch_id(s, fsw, a);
				int id_f = switch_id(s, lo_w, a);
				swleaving += x[id_fsw] * turn_SW;
				fhr += x[id_f];
			}
		}
		turn_hr = swleaving / fhr;
		for (int s = 0; s < S; s++) {
			for (int a = 0; a < A; a++) {
				int id_fsw = switch_id(s, fsw, a);
				int id_f = switch_id(s, lo_w, a);
				turnover_ou[lo_w][s][a] = -x[id_f] * (turn_hr);
				turnover_ou[fsw][s][a] = -x[id_fsw] * turn_SW;

				turnover_in[lo_w][s][a] = x[id_fsw] * turn_SW;
				turnover_in[fsw][s][a] = x[id_f] * turn_hr;
			}
		}
		for (int k = 0; k < R; k++) {
			for (int s = 0; s < S; s++) {
				for (int a = 0; a < A; a++) {
					turnover[k][s][a] = turnover_in[k][s][a] + turnover_ou[k][s][a];
				}
			}
		}

		//.................................................................................................
		// >>>>>>>>>>>>>>>>>>>>>>>> ODE's <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
		//................................................................................................. 

		// Outcome temporal storage

		double out[noutcomes] = { 0.0 };
		int g = nstages;// switch_id(S, R, A);

		for (int r = 0; r < R; r++) {
			for (int a = 0; a < A; a++) {







				//Define Indexes for all S by specific a and r
				int iS = Sus * R * A + r * A + a;
				int iE = Ac * R * A + r * A + a;
				int iU1 = U_500p * R * A + r * A + a;
				int iU2 = U_500 * R * A + r * A + a;
				int iU3 = U_350 * R * A + r * A + a;
				int iU4 = U_200 * R * A + r * A + a;
				int iD1 = D_500p * R * A + r * A + a;
				int iD2 = D_500 * R * A + r * A + a;
				int iD3 = D_350 * R * A + r * A + a;
				int iD4 = D_200 * R * A + r * A + a;
				int iA1 = A_500p * R * A + r * A + a;
				int iA2 = A_500 * R * A + r * A + a;
				int iA3 = A_350 * R * A + r * A + a;
				int iA4 = A_200 * R * A + r * A + a;
				int iS1 = S_500p * R * A + r * A + a;
				int iS2 = S_500 * R * A + r * A + a;
				int iS3 = S_350 * R * A + r * A + a;
				int iS4 = S_200 * R * A + r * A + a;
				int iDi1 = Di_500p * R * A + r * A + a;
				int iDi2 = Di_500 * R * A + r * A + a;
				int iDi3 = Di_350 * R * A + r * A + a;
				int iDi4 = Di_200 * R * A + r * A + a;

				// Susceptibles (other indexes on s and d for susceptibles wont be used)

				if (a == 0) {
					dxdt[iS] = Recruit[r] - x[iS] * (FOI_grand[r][a] + mu[r][a]) + turnover[r][Sus][a];
				}
				else
				{
					dxdt[iS] = turnover[r][Sus][a] - x[iS] * (FOI_grand[r][a] + mu[r][a]);
				}



				double test_adj = 1;// test*Tratios[r][a];

				/*if ((fmod(t, 1) == 0 || t == dt) && t != 0) {
					cout << t << " " << x[iD1] << " " << x[iD2] << " " << x[iD3] << " " << x[iD4] << endl;
				}
				*/

				// Acute/Early HIv Infection 
				dxdt[iE] = x[iS] * FOI_grand[r][a] - x[iE] * (proghorz[0] + mu[r][a]) + turnover[r][Ac][a];
				//Undetected stages 
				dxdt[iU1] = x[iE] * (proghorz[0] * postAcute[0]) - x[iU1] * (test_adj*hivT_scale + proghorz[1] + mu[r][a] + mortHIV[0]) + turnover[r][U_500p][a];
				dxdt[iU2] = x[iE] * (proghorz[0] * postAcute[1]) + x[iU1] * proghorz[1] - x[iU2] * (test_adj*hivT_scale + proghorz[2] + mu[r][a] + mortHIV[1]) + turnover[r][U_500][a];
				dxdt[iU3] = x[iE] * (proghorz[0] * postAcute[2]) + x[iU2] * proghorz[2] - x[iU3] * (test *test_adj + proghorz[3] + mu[r][a] + mortHIV[2]) + turnover[r][U_350][a];
				dxdt[iU4] = x[iE] * (proghorz[0] * postAcute[3]) + x[iU3] * proghorz[3] - x[iU4] * (AIDStestRR*test *test_adj + mortHIV[3] + mu[r][a]) + turnover[r][U_200][a];
				//Detected stages 
				dxdt[iD1] = x[iU1] * test_adj*hivT_scale * (1 - p_art) - x[iD1] * (mortHIV[0] + proghorz[1] + treatment[0] + mu[r][a]) + turnover[r][D_500p][a];
				dxdt[iD2] = x[iU2] * test_adj*hivT_scale * (1 - p_art) + x[iD1] * proghorz[1] - x[iD2] * (mortHIV[1] + proghorz[2] + treatment[1] + mu[r][a]) + turnover[r][D_500][a];
				dxdt[iD3] = x[iU3] * test_adj * test * (1 - p_art) + x[iD2] * proghorz[2] - x[iD3] * (mortHIV[2] + proghorz[3] + treatment[2] + mu[r][a]) + turnover[r][D_350][a];
				dxdt[iD4] = x[iU4] * AIDStestRR * test * test_adj * (1 - p_art) + x[iD3] * proghorz[3] - x[iD4] * (mortHIV[3] + treatment[3] + mu[r][a]) + turnover[r][D_200][a];
				//ART stages 
				dxdt[iA1] = x[iU1] * test_adj*hivT_scale * p_art + x[iD1] * treatment[0] + x[iDi1] * back_art[0] - x[iA1] * (supp + mortSupp[0] + mu[r][a]) + turnover[r][A_500p][a];
				dxdt[iA2] = x[iU2] * test_adj*hivT_scale * p_art + x[iD2] * treatment[1] + x[iDi2] * back_art[1] - x[iA2] * (supp + mortSupp[1] + mu[r][a]) + turnover[r][A_500][a];
				dxdt[iA3] = x[iU3] * test *test_adj * p_art + x[iD3] * treatment[2] + x[iDi3] * back_art[2] - x[iA3] * (supp + mortSupp[2] + mu[r][a]) + turnover[r][A_350][a];
				dxdt[iA4] = x[iU4] * AIDStestRR * test *test_adj * p_art + x[iD4] * treatment[3] + x[iDi4] * back_art[3] - x[iA4] * (supp + mortSupp[3] + mu[r][a]) + turnover[r][A_200][a];
				//Suppressed stages 
				dxdt[iS1] = x[iA1] * supp - x[iS1] * (dis + mortSupp[0] + mu[r][a]) + turnover[r][S_500p][a];
				dxdt[iS2] = x[iA2] * supp - x[iS2] * (dis + mortSupp[1] + mu[r][a]) + turnover[r][S_500][a];
				dxdt[iS3] = x[iA3] * supp - x[iS3] * (dis + mortSupp[2] + mu[r][a]) + turnover[r][S_350][a];
				dxdt[iS4] = x[iA4] * supp - x[iS4] * (dis + mortSupp[3] + mu[r][a]) + turnover[r][S_200][a];
				//Disengaged stages 
				dxdt[iDi1] = x[iS1] * dis - x[iDi1] * (back_art[0] + proghorz[1] + mortHIV[0] + mu[r][a]) + turnover[r][Di_500p][a];
				dxdt[iDi2] = x[iS2] * dis + x[iDi1] * proghorz[1] - x[iDi2] * (back_art[1] + proghorz[2] + mortHIV[1] + mu[r][a]) + turnover[r][Di_500][a];
				dxdt[iDi3] = x[iS3] * dis + x[iDi2] * proghorz[2] - x[iDi3] * (back_art[2] + proghorz[3] + mortHIV[2] + mu[r][a]) + turnover[r][Di_350][a];
				dxdt[iDi4] = x[iS4] * dis + x[iDi3] * proghorz[3] - x[iDi4] * (back_art[3] + mortHIV[3] + mu[r][a]) + turnover[r][Di_200][a];

				// Temp cumultive outcomes
				out[prev - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
				out[inc - g] += x[iS] * FOI_grand[r][a];
				out[sus - g] += x[iS];
				out[undet - g] += x[iU1] + x[iU2] + x[iU3] + x[iU4];
				out[detec - g] += x[iD1] + x[iD2] + x[iD3] + x[iD4];
				out[onart - g] += x[iA1] + x[iA2] + x[iA3] + x[iA4];
				out[suppr - g] += x[iS1] + x[iS2] + x[iS3] + x[iS4];
				out[diseng - g] += x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
				out[mortundet - g] += x[iU1] * mortHIV[0] + x[iU2] * mortHIV[1] + x[iU3] * mortHIV[2] + x[iU4] * mortHIV[3];
				out[mortdetec - g] += x[iD1] * mortHIV[0] + x[iD2] * mortHIV[1] + x[iD3] * mortHIV[2] + x[iD4] * mortHIV[3];
				out[mortonart - g] += x[iA1] * mortSupp[0] + x[iA2] * mortSupp[2] + x[iA3] * mortSupp[2] + x[iA4] * mortSupp[3];
				out[mortsuppr - g] += x[iS1] * mortSupp[0] + x[iS2] * mortSupp[1] + x[iS3] * mortSupp[2] + x[iS4] * mortSupp[3];
				out[mortdiseng - g] += x[iDi1] * mortHIV[0] + x[iDi2] * mortHIV[1] + x[iDi3] * mortHIV[2] + x[iDi4] * mortHIV[3];
				out[mortaHIV - g] += x[iU1] * mortHIV[0] + x[iU2] * mortHIV[1] + x[iU3] * mortHIV[2] + x[iU4] * mortHIV[3] + x[iD1] * mortHIV[0] + x[iD2] * mortHIV[1] + x[iD3] * mortHIV[2] + x[iD4] * mortHIV[3] +
					x[iA1] * mortHIV[0] + x[iA2] * mortHIV[2] + x[iA3] * mortHIV[2] + x[iA4] * mortHIV[3] + x[iS1] * mortSupp[0] + x[iS2] * mortSupp[1] + x[iS3] * mortSupp[2] + x[iS4] * mortSupp[3] +
					x[iDi1] * mortHIV[0] + x[iDi2] * mortHIV[1] + x[iDi3] * mortHIV[2] + x[iDi4] * mortHIV[3];
				out[mortaAIDS - g] += (x[iD4] * mortHIV[3] + x[iA4] * mortSupp[3] + x[iS4] * mortSupp[3] + x[iDi4] * mortHIV[3]) + x[iU4] * mortHIV[3] * undeR;
				out[newdetec - g] += (x[iU1] + x[iU2])*test_adj*hivT_scale + (x[iU4] * AIDStestRR * test_adj) + x[iU3] * test_adj;
				out[newart - g] += (x[iU1] + x[iU2])*test_adj*hivT_scale*p_art + x[iU3] * test_adj*p_art + x[iU4] * test_adj *p_art*AIDStestRR + x[iD1] * treatment[0] + x[iD2] * treatment[1] + x[iD3] * treatment[2] + x[iD4] * treatment[3];
				out[dxdeath - g] += (x[iU1] * mortHIV[0] + x[iU2] * mortHIV[1] + x[iU3] * mortHIV[2] + x[iU4] * mortHIV[3]) * undeR;
				out[dxaids - g] += x[iU4] * AIDStestRR * test*test_adj + x[iU3] * test*test_adj;
				out[dxhiv - g] += (x[iU1] + x[iU2]) * test_adj*hivT_scale;
				out[plha - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];

				if (r == lo_w) {
					out[prevlo_w - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					out[inclo_w - g] += x[iS] * FOI_grand[r][a];
					out[suslo_w - g] += x[iS];
					out[poplo_w - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
				}
				if (r == fsw) {
					out[prevfsw - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					out[incfsw - g] += x[iS] * FOI_grand[r][a];
					out[susfsw - g] += x[iS];
					out[popfsw - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
				}
				if (r == lo_m) {
					out[prevlo_m - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					out[inclo_m - g] += x[iS] * FOI_grand[r][a];
					out[suslo_m - g] += x[iS];
					out[poplo_m - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
				}
				if (r == cfsw) {
					out[prevcfsw - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					out[inccfsw - g] += x[iS] * FOI_grand[r][a];
					out[suscfsw - g] += x[iS];
					out[popcfsw - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
				}
				if (r == bi) {
					out[prevbi - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					out[incbi - g] += x[iS] * FOI_grand[r][a];
					out[susbi - g] += x[iS];
					out[popbi - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
				}
				if (r == msm) {
					out[prevmsm - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					out[incmsm - g] += x[iS] * FOI_grand[r][a];
					out[susmsm - g] += x[iS];
					out[popmsm - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
				}
				if (a == 0) {
					out[notaids_a0 - g] += x[iU4] * test* AIDStestRR * test_adj + x[iU3] * test* test_adj;
					out[nothiv_a0 - g] += (x[iU1] + x[iU2]) * test_adj*hivT_scale;
					out[notdead_a0 - g] += x[iU4] * mortHIV[3] * undeR;
					if (r == lo_w || r == fsw) {
						out[prevpreg15 - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
						out[poppreg15 - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					}
					if (r == nr_m || r == lo_m || r == cfsw || r == bi || r == msm) {
						out[prevM15 - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
						out[popM15 - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					}
				}
				if (a == 1) {
					out[notaids_a1 - g] += x[iU4] * AIDStestRR * test * test_adj + x[iU3] * test * test_adj;
					out[nothiv_a1 - g] += (x[iU1] + x[iU2]) * test_adj*hivT_scale;
					out[notdead_a1 - g] += x[iU4] * mortHIV[3] * undeR;
					if (r == lo_w || r == fsw) {
						out[prevpreg34 - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
						out[poppreg34 - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					}
					if (r == nr_m || r == lo_m || r == cfsw || r == bi || r == msm) {
						out[prevM34 - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
						out[popM34 - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					}
				}
				if (a == 2) {
					out[notaids_a2 - g] += x[iU4] * AIDStestRR * test_adj * test + x[iU3] * test* test_adj;
					out[nothiv_a2 - g] += (x[iU1] + x[iU2]) * test_adj*hivT_scale;
					out[notdead_a2 - g] += x[iU4] * mortHIV[3] * undeR;
				}
				if (a == 3) {
					out[notaids_a3 - g] += x[iU4] * AIDStestRR * test * test_adj + x[iU3] * test* test_adj;
					out[nothiv_a3 - g] += (x[iU1] + x[iU2]) * test_adj*hivT_scale;
					out[notdead_a3 - g] += x[iU4] * mortHIV[3] * undeR;
				}
				if (a < 2) {
					if (r == lo_w || r == fsw) {
						out[prevpreg49 - g] += x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
						out[poppreg49 - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					}
					/*if (r == lo_w || r == fsw || r==nr_w) {
						out[poppreg49 - g] += x[iS] + x[iE] + x[iU1] + x[iU2] + x[iU3] + x[iU4] + x[iD1] + x[iD2] + x[iD3] + x[iD4] + x[iA1] + x[iA2] + x[iA3] + x[iA4] + x[iS1] + x[iS2] + x[iS3] + x[iS4] + x[iDi1] + x[iDi2] + x[iDi3] + x[iDi4];
					}*/


				}

			}//a loop
		}//r loop


				///////////////////////////// OUTCOMES ///////////////////////////////////////////////
				//// Fill rest of structure with output
		out[pop - g] = totPop;



		for (int jj = g; jj < S_all; jj++) {

			dxdt[jj] = out[jj - g];

		}

		if (t > timeZero) {
			dxdt[lifeyears] = x[lifeyears] + x[pop];
		}

		delete[] HIVratio;
		delete[] Dis500p;
		delete[] Dis500;
		delete[] Dis350;
		delete[] Dis200;
		delete[] Back500p;
		delete[] Back500;
		delete[] Back350;
		delete[] Back200;
		delete[] ART350p;
		delete[] ART350;
		delete[] ART200;
		delete[] Supp500p;
		delete[] Supp500;
		delete[] Supp350;
		delete[] Supp200;
		delete[] proptoART;

	}// End of operator function
}; // end of ODE's class'
//--------------------------------------------------------------------------------------------------------------------------------------------


//--------------------Function for output collection at every unit time (yr)-------------------------------------------------------------------

//class Write_results
struct Write_results
{
	Results * res;
	int time_0;
	int time_end;
	size_t m_steps;
	size_t max_steps;
	//
	//public:
	//Write_results(Results * Re, int time_p) : res(Re), time_end(time_p) {}
	Write_results(Results * Re, int time_i, int time_p, size_t maxS) : res(Re), time_0(time_i), time_end(time_p), m_steps(0), max_steps(maxS) {}
	// write results function
	void operator()(const state_type &x, const double t)
	{

		m_steps++;
		if (m_steps > max_steps) throw runtime_error("Too much steps");


		//// Pass dxdt contents for 47 stages

		if (fmod(t, 1.0) == 0.0) {// collect results every year
			int ii = int(t) - time_0;
			for (int j = 0; j < S_all; j++) {
				res->state[j][ii] = x[j];
			};

			int Dursim = 1 + (time_end - time_0);

			if (t == time_end) {
				// get the differential of outcome (cumulative) stages
				vector<double> prev_all, prev_m, prev_sw, prev_f, prev_cfsw, prev_bi, prev_xmsm, prev_msm, prev_preg15, prev_preg34, prev_preg49, prev_m15, prev_m34,
					inc_all, inc_m, inc_sw, inc_f, inc_cfsw, inc_bi, inc_xmsm, inc_msm, sus_all, sus_m, sus_sw, sus_f, sus_cfsw, sus_bi, sus_xmsm, sus_msm,
					undetected, detected, onART, suppresed, disengaged, mort_undet, mort_detec, mort_onart, mort_suppr, mort_diseng, mort_HIV, mort_AIDS, new_detec, new_art, dx_death, dx_aids, dx_hiv,
					lifeyrs, PLHA, pop_all, pop_m, pop_sw, pop_f, pop_cfsw, pop_bi, pop_xmsm, pop_msm, pop_preg15, pop_preg34, pop_preg49, pop_m15, pop_m34,
					notA0, notA1, notA2, notA3, notH0, notH1, notH2, notH3, notD0, notD1, notD2, notD3;
				for (int ii = 0; ii < Dursim; ii++) {
					prev_all.push_back(res->state[prev][ii]);
					prev_f.push_back(res->state[prevlo_w][ii]);
					prev_sw.push_back(res->state[prevfsw][ii]);
					prev_m.push_back(res->state[prevlo_m][ii]);
					prev_cfsw.push_back(res->state[prevcfsw][ii]);
					prev_bi.push_back(res->state[prevbi][ii]);
					prev_xmsm.push_back(res->state[prevmsm][ii]);
					prev_msm.push_back((res->state[prevmsm][ii]) + (res->state[prevbi][ii]));///
					prev_preg15.push_back(res->state[prevpreg15][ii]);
					prev_preg34.push_back(res->state[prevpreg34][ii]);
					prev_preg49.push_back(res->state[prevpreg49][ii]  + res->state[prevpreg34][ii]+ res->state[prevpreg15][ii]);
					prev_m15.push_back(res->state[prevM15][ii]);
					prev_m34.push_back(res->state[prevM34][ii]);
					inc_all.push_back(res->state[inc][ii]);
					inc_f.push_back(res->state[inclo_w][ii]);
					inc_sw.push_back(res->state[incfsw][ii]);
					inc_m.push_back(res->state[inclo_m][ii]);
					inc_cfsw.push_back(res->state[inccfsw][ii]);
					inc_bi.push_back(res->state[incbi][ii]);
					inc_xmsm.push_back(res->state[incmsm][ii]);
					inc_msm.push_back((res->state[incmsm][ii]) + (res->state[incbi][ii]));///
					sus_all.push_back(res->state[sus][ii]);
					sus_f.push_back(res->state[suslo_w][ii]);
					sus_sw.push_back(res->state[susfsw][ii]);
					sus_m.push_back(res->state[suslo_m][ii]);
					sus_cfsw.push_back(res->state[suscfsw][ii]);
					sus_bi.push_back(res->state[susbi][ii]);
					sus_xmsm.push_back(res->state[susmsm][ii]);
					sus_msm.push_back((res->state[susmsm][ii]) + (res->state[susbi][ii]));///
					undetected.push_back(res->state[undet][ii]);
					detected.push_back(res->state[detec][ii]);
					onART.push_back(res->state[onart][ii]);
					suppresed.push_back(res->state[suppr][ii]);
					disengaged.push_back(res->state[diseng][ii]);
					mort_undet.push_back(res->state[mortundet][ii]);
					mort_detec.push_back(res->state[mortdetec][ii]);
					mort_onart.push_back(res->state[mortonart][ii]);
					mort_suppr.push_back(res->state[mortsuppr][ii]);
					mort_diseng.push_back(res->state[mortdiseng][ii]);
					mort_HIV.push_back(res->state[mortaHIV][ii]);
					mort_AIDS.push_back(res->state[mortaAIDS][ii]);
					new_detec.push_back(res->state[newdetec][ii]);
					new_art.push_back(res->state[newart][ii]);
					dx_death.push_back(res->state[dxdeath][ii]);
					dx_aids.push_back(res->state[dxaids][ii]);
					dx_hiv.push_back(res->state[dxhiv][ii]);
					lifeyrs.push_back(res->state[lifeyears][ii]);
					PLHA.push_back(res->state[plha][ii]);
					pop_all.push_back(res->state[pop][ii]);
					pop_f.push_back(res->state[poplo_w][ii]);
					pop_sw.push_back(res->state[popfsw][ii]);
					pop_m.push_back(res->state[poplo_m][ii]);
					pop_cfsw.push_back(res->state[popcfsw][ii]);
					pop_bi.push_back(res->state[popbi][ii]);
					pop_xmsm.push_back(res->state[popmsm][ii]);
					pop_msm.push_back((res->state[popmsm][ii]) + (res->state[popbi][ii]));///
					pop_preg15.push_back(res->state[poppreg15][ii]);
					pop_preg34.push_back(res->state[poppreg34][ii]);
					pop_preg49.push_back(res->state[poppreg49][ii] + res->state[poppreg34][ii]+ res->state[poppreg15][ii]);
					pop_m15.push_back(res->state[popM15][ii]);
					pop_m34.push_back(res->state[popM34][ii]);
					notA0.push_back(res->state[notaids_a0][ii]);
					notA1.push_back(res->state[notaids_a1][ii]);
					notA2.push_back(res->state[notaids_a2][ii]);
					notA3.push_back(res->state[notaids_a3][ii]);
					notH0.push_back(res->state[nothiv_a0][ii]);
					notH1.push_back(res->state[nothiv_a1][ii]);
					notH2.push_back(res->state[nothiv_a2][ii]);
					notH3.push_back(res->state[nothiv_a3][ii]);
					notD0.push_back(res->state[notdead_a0][ii]);
					notD1.push_back(res->state[notdead_a1][ii]);
					notD2.push_back(res->state[notdead_a2][ii]);
					notD3.push_back(res->state[notdead_a3][ii]);

				}

				prev_all = diff(prev_all);
				prev_f = diff(prev_f);
				prev_sw = diff(prev_sw);
				prev_m = diff(prev_m);
				prev_cfsw = diff(prev_cfsw);
				prev_bi = diff(prev_bi);
				prev_xmsm = diff(prev_xmsm);
				prev_msm = diff(prev_msm);
				prev_preg15 = diff(prev_preg15);
				prev_preg34 = diff(prev_preg34);
				prev_preg49 = diff(prev_preg49);
				prev_m15 = diff(prev_m15);
				prev_m34 = diff(prev_m34);
				inc_all = diff(inc_all);
				inc_f = diff(inc_f);
				inc_sw = diff(inc_sw);
				inc_m = diff(inc_m);
				inc_cfsw = diff(inc_cfsw);
				inc_bi = diff(inc_bi);
				inc_xmsm = diff(inc_xmsm);
				inc_msm = diff(inc_msm);
				sus_all = diff(sus_all);
				sus_f = diff(sus_f);
				sus_sw = diff(sus_sw);
				sus_m = diff(sus_m);
				sus_cfsw = diff(sus_cfsw);
				sus_bi = diff(sus_bi);
				sus_xmsm = diff(sus_xmsm);
				sus_msm = diff(sus_msm);
				undetected = diff(undetected);
				detected = diff(detected);
				onART = diff(onART);
				suppresed = diff(suppresed);
				disengaged = diff(disengaged);
				mort_undet = diff(mort_undet);
				mort_detec = diff(mort_detec);
				mort_onart = diff(mort_onart);
				mort_suppr = diff(mort_suppr);
				mort_diseng = diff(mort_diseng);
				mort_HIV = diff(mort_HIV);
				mort_AIDS = diff(mort_AIDS);
				new_detec = diff(new_detec);
				new_art = diff(new_art);
				dx_death = diff(dx_death);
				dx_aids = diff(dx_aids);
				dx_hiv = diff(dx_hiv);
				lifeyrs = diff(lifeyrs);
				PLHA = diff(PLHA);
				pop_all = diff(pop_all);
				pop_f = diff(pop_f);
				pop_sw = diff(pop_sw);
				pop_m = diff(pop_m);
				pop_cfsw = diff(pop_cfsw);
				pop_bi = diff(pop_bi);
				pop_xmsm = diff(pop_xmsm);
				pop_msm = diff(pop_msm);
				pop_preg15 = diff(pop_preg15);
				pop_preg34 = diff(pop_preg34);
				pop_preg49 = diff(pop_preg49);
				pop_m15 = diff(pop_m15);
				pop_m34 = diff(pop_m34);
				notA0 = diff(notA0);
				notA1 = diff(notA1);
				notA2 = diff(notA2);
				notA3 = diff(notA3);
				notD0 = diff(notD0);
				notD1 = diff(notD1);
				notD2 = diff(notD2);
				notD3 = diff(notH3);
				notH0 = diff(notH0);
				notH1 = diff(notH1);
				notH2 = diff(notH2);
				notH3 = diff(notH3);
				for (int h = 0; h < Dursim; h++) {
					res->state[prev][h] = prev_all[h] / pop_all[h];  //cout << prev_all[h] << endl;
					res->state[prevlo_w][h] = prev_f[h] / pop_f[h];
					res->state[prevfsw][h] = prev_sw[h] / pop_sw[h];
					res->state[prevlo_m][h] = prev_m[h] / pop_m[h];
					res->state[prevcfsw][h] = prev_cfsw[h] / pop_cfsw[h];
					res->state[prevbi][h] = prev_bi[h] / pop_bi[h];
					res->state[prevmsm][h] = prev_xmsm[h] / pop_xmsm[h];
					res->state[prevmsmall][h] = prev_msm[h] / pop_msm[h];
					res->state[prevpreg15][h] = prev_preg15[h] / pop_preg15[h]; //cout <<"prevpreg15 "<< prev_preg15[h] / pop_preg15[h] << endl;
					res->state[prevpreg34][h] = prev_preg34[h] / pop_preg34[h];
					res->state[prevpreg49][h] = prev_preg49[h] / pop_preg49[h];
					res->state[prevM15][h] = prev_m15[h] / pop_m15[h];
					res->state[prevM34][h] = prev_m34[h] / pop_m34[h];
					res->state[inc][h] = inc_all[h] / sus_all[h];
					res->state[inclo_w][h] = inc_f[h] / sus_f[h];
					res->state[incfsw][h] = inc_sw[h] / sus_sw[h];
					res->state[inclo_m][h] = inc_m[h] / sus_m[h];
					res->state[inccfsw][h] = inc_cfsw[h] / sus_cfsw[h];
					res->state[incbi][h] = inc_bi[h] / sus_bi[h];
					res->state[incmsm][h] = inc_xmsm[h] / sus_xmsm[h];
					res->state[incmsmall][h] = inc_msm[h] / sus_msm[h];
					res->state[newinf][h] = inc_all[h];
					res->state[newinflo_w][h] = inc_f[h];
					res->state[newinffsw][h] = inc_sw[h];
					res->state[newinflo_m][h] = inc_m[h];
					res->state[newinfcfsw][h] = inc_cfsw[h];
					res->state[newinfbi][h] = inc_bi[h];
					res->state[newinfmsm][h] = inc_xmsm[h];
					res->state[newinfmsmall][h] = inc_msm[h];
					res->state[undet][h] = undetected[h];
					res->state[detec][h] = detected[h];
					res->state[onart][h] = onART[h];  // cout <<onART[h] << endl;
					res->state[suppr][h] = suppresed[h];
					res->state[diseng][h] = disengaged[h];
					res->state[mortundet][h] = mort_undet[h];
					res->state[mortdetec][h] = mort_detec[h];
					res->state[mortonart][h] = mort_onart[h];
					res->state[mortsuppr][h] = mort_suppr[h];
					res->state[mortdiseng][h] = mort_diseng[h];
					res->state[mortaHIV][h] = mort_HIV[h];
					res->state[mortaAIDS][h] = mort_AIDS[h];
					res->state[newdetec][h] = new_detec[h];
					res->state[newart][h] = new_art[h];
					res->state[dxdeath][h] = dx_death[h];
					res->state[dxaids][h] = dx_aids[h];
					res->state[dxhiv][h] = dx_hiv[h];
					res->state[lifeyears][h] = lifeyrs[h];
					res->state[plha][h] = PLHA[h];
					res->state[pop][h] = pop_all[h];
					res->state[notaids_a0][h] = notA0[h];
					res->state[notaids_a1][h] = notA1[h];
					res->state[notaids_a2][h] = notA2[h];
					res->state[notaids_a3][h] = notA3[h];
					res->state[nothiv_a0][h] = notH0[h];
					res->state[nothiv_a1][h] = notH1[h];
					res->state[nothiv_a2][h] = notH2[h];
					res->state[nothiv_a3][h] = notH3[h];
					res->state[notdead_a0][h] = notD0[h];
					res->state[notdead_a1][h] = notD1[h];
					res->state[notdead_a2][h] = notD2[h];
					res->state[notdead_a3][h] = notD3[h];

				}

				//fix nan at t0
				res->state[prev][0] = 0;
				res->state[prevlo_w][0] = 0;
				res->state[prevfsw][0] = 0;
				res->state[prevlo_m][0] = 0;
				res->state[prevcfsw][0] = 0;
				res->state[prevbi][0] = 0;
				res->state[prevmsm][0] = 0;
				res->state[prevmsmall][0] = 0;
				res->state[prevpreg15][0] = 0;
				res->state[prevpreg34][0] = 0;
				res->state[prevpreg49][0] = 0;
				res->state[prevM15][0] = 0;
				res->state[prevM34][0] = 0;
				res->state[inc][0] = 0;
				res->state[inclo_w][0] = 0;
				res->state[incfsw][0] = 0;
				res->state[inclo_m][0] = 0;
				res->state[inccfsw][0] = 0;
				res->state[incbi][0] = 0;
				res->state[incmsm][0] = 0;
				res->state[incmsmall][0] = 0;

			}
		}

	}
};

//-------------------------------------------------------------------------------------------------------------------------------------------------



//--------------------------------------// Model Class wrapper (runs model )--------------------------------------------------------------------------

Model::Model()
{
	//ctor
}

Model::~Model()
{
	//dtor
}

Results * Model::Run(double * Importedparams, int sizep, int time0, int time_end, std::vector<double> interv) {  // Wraping Results function

	double seed = Importedparams[0]; //0.0005;//
	double FSW_frac = Importedparams[1]; //0.0005;//
	double CFSW_frac = Importedparams[2]; //0.0005;//
	double MSM_frac = Importedparams[3]; //0.0005;//

	state_type x;// Create mutidim array and zero it
	for (int s = 0; s < S_all; s++) {
		x[s] = 0.0;
	}

	// Distribution of Risk categories

	double Q[R] = { 0.0 };
	Q[nr_w] = NoRW_frac * f_prop;//proportion of no risk females
	Q[lo_w] = (1 - (FSW_frac + NoRW_frac))*f_prop;// proportion Non-FSW Women
	Q[fsw] = FSW_frac * f_prop;// proportion of female sex Workers

	Q[nr_m] = NoRM_frac * m_prop;// proportion of no risk males
	Q[lo_m] = (1 - (NoRM_frac + CFSW_frac + MSM_frac))*m_prop; // proportion of heterosexual men
	Q[cfsw] = CFSW_frac * m_prop;// proportion of clients of female sex Workers
	Q[bi] = Bi_frac*MSM_frac * m_prop; // proportion of low risk MSM
	Q[msm] = (1 - Bi_frac)*MSM_frac * m_prop; // proportion of high risk MSM



	// Initial Pop
	double ninit = initPop;

	// Seed in risk categories and set initial pop disrtibution
	for (int a = 0; a < A; a++) {


		//x[Ac * R * A + 1 * A + a] = ninit*(Q[1] * ageDist[a])*seed;//lo
		x[Ac * R * A + 2 * A + a] = ninit*(Q[2] * ageDist[a])*seed;//fsw
		//x[Ac * R * A + 4 * A + a] = ninit*(Q[4] * ageDist[a])*seed;//lo
		//x[Ac * R * A + 5 * A + a] = ninit*(Q[5] * ageDist[a])*seed;//cfsw
		//x[Ac * R * A + 6 * A + a] = ninit*(Q[6] * ageDist[a])*seed;//Bi
		x[Ac * R * A + 7 * A + a] = ninit*(Q[7] * ageDist[a])*seed;//MSM



		x[Sus * R * A + 0 * A + a] = ninit*(Q[0] * ageDist[a]) - x[Ac * R * A + 0 * A + a];
		x[Sus * R * A + 1 * A + a] = ninit*(Q[1] * ageDist[a]) - x[Ac * R * A + 1 * A + a];
		x[Sus * R * A + 2 * A + a] = ninit*(Q[2] * ageDist[a]) - x[Ac * R * A + 2 * A + a];
		x[Sus * R * A + 3 * A + a] = ninit*(Q[3] * ageDist[a]) - x[Ac * R * A + 3 * A + a];
		x[Sus * R * A + 4 * A + a] = ninit*(Q[4] * ageDist[a]) - x[Ac * R * A + 4 * A + a];
		x[Sus * R * A + 5 * A + a] = ninit*(Q[5] * ageDist[a]) - x[Ac * R * A + 5 * A + a];
		x[Sus * R * A + 6 * A + a] = ninit*(Q[6] * ageDist[a]) - x[Ac * R * A + 6 * A + a];
		x[Sus * R * A + 7 * A + a] = ninit*(Q[7] * ageDist[a]) - x[Ac * R * A + 7 * A + a];

	}

	//typedef runge_kutta_dopri5<state_type> rkdp;


	typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_dopri5< state_type > > dopri_stepper_type;
	typedef boost::numeric::odeint::dense_output_runge_kutta< dopri_stepper_type > dense_stepper_type;


	int length = 1 + (time_end - timeZero);
	Results * ModelResults = new Results(length);    //Create structure that contains Results

	//    int sizeiv = sizeof(interv)/sizeof(interv[0]);
	input pars;    // Create parameter vector
	for (int p = 0; p < sizep; p++) {
		pars.push_back(Importedparams[p]);
	}

	for (int pp = 0; pp < interv.size(); pp++) {
		pars.push_back(interv[pp]);
	}

	odes brm(pars); // create structure with input params and ODEs
	size_t steps = 0;
	size_t max_steps = (length)* 100;
	Write_results results_fn(ModelResults, time0, time_end, max_steps); // Create object that will read and write results
	try {
		//  steps=integrate_const(make_controlled(1E-6, 1E-6, rkdp()), brm, x, 0.0, double(time_end), dt, results_fn);// rungeKutta dormadprince

		steps = integrate_const(dense_stepper_type(), brm, x, double(time0), double(time_end), dt, results_fn);
	}
	catch (...) { steps = max_steps; }

	return ModelResults;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------
