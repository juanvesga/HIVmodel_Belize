//
//  dataheader.h
//  MCMCsoutheast
//
//  Created by juan fernando vesga on 21/03/2016.
//  Copyright © 2016 juan fernando vesga. All rights reserved.
//

#ifndef dataheader_h
#define dataheader_h
#include <iostream>
#include <vector>
using namespace std;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::Calibration Data::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



// ------------------------------------------------Parameter Bounds for sampling
double bounds[24][2] = {

	{ 1e-6,	1e-3 },//	seed
	{ 0.002,	0.01 },//  FSWfrac
	{ 0.05,	    0.12 },//	CFSWfrac
	{ 0.01,	    0.1 },//	MSMfrac"
	{ 0,        1 },    //	epsi
	{ 0.00013,	0.00141 },//betaFtoM"
	{ 0.0006,	0.00109 },//betaMtoF"
	{2         , 20},// { 0.0006,	0.025 },//	BetaMSM	{ 0.002,	0.025 },//	BetaMSM
	{ 0.1,	    0.6 },// detectionrpryr
	{ 2000,	    2017 },//		test_peak
	{ 0.01,	    0.5 },//	ART_coef"
    {2   ,       20},// c_msm = prm[11];
    {1   ,       10 },//c_bi = prm[12];
    {2    ,      20},//c_cfsw = prm[13];
    {120 ,      300 },//c_fsw = prm[14];
	{ 0 ,     4 },//c_m = prm[14];
	{ 0 ,     4 },//c_f = prm[14];
	{ 2 ,     20 },//actsMSM = prm[14];
	{1 ,      10 },//condom_x
    {0.5 ,      1 },//propBiacts
    {0.5 ,      1.5 },//condymsm 
	{ 0.5 ,      1.5 },//condysw 
	{ 0.5 ,      1.5 },//condyhet 
    {4.7 ,     18.8}//Acute coefficient



//{ 0.00001,	0.0001 },//	seed
//{ 0.002,	0.005 },//  FSWfrac
//{ 0.05,	    0.1 },//	CFSWfrac
//{ 0.01,	    0.08 },//	MSMfrac"
//{ 0,        1 },    //	epsi
//{ 0.00013,	0.00141 },//betaFtoM"
//{ 0.0006,	0.00109 },//betaMtoF"
//{ 0.002,	0.025 },//	BetaMSM
//{ 0.1,	    1 },// detectionrpryr
//{ 1,	    10 },//		sigma_test
//{ 0.01,	    0.5 },//	ART_coef"
//{ 5   ,        15 },// c_msm = prm[11];
//{ 5   ,       15 },//c_bi = prm[12];
//{ 3    ,       10 },//c_cfsw = prm[13];
//{ 120 ,      250 }//c_fsw = prm[14];

};

// >>>>>>>>>>>>>>>>>>>>>>>> ALL MODEL DATA
		vector<int> gen_years{ 2007, 2014 };
		vector<int> msm_years{ 2014 };
		vector<int> fsw_years{ 2012 };
		vector<int> preg_years{ 2006,2007,2008,2009,2010 };
		vector<int> hivcases_years{ 1986,	1987,	1988,	1989,	1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010 };
		vector<int> onART_years = { 2007,2008,2009,2010 };
		vector<int> mort_years = { 2008,	2009,	2010,	2011	,2012,	2013,	2014,	2015,	2016 };



		double prevgen_data[2] = { 0.021, 0.014 };
		double prevmsm_data= 0.1385;
		double prevfsw_data = 0.0091;
		double prevpreg_data[5] = { 0.0097,	0.0098,	0.0099,	0.0099,	0.0086 };
		double hivcases_data[25]{ 1,	13,	19,	21,	34,	44,	75,	81,	60,	72,	93,	99,	184,	241,	226,	330,	431,	447,	457,	434,	443,	450,	425,	365,	244};
		double onART_data[4] = { 558, 630,	855, 1053 };
		double mort_data[9] = { 62,	106,	114,	99,	109,	110,	100,	113,	104 };
		
		// Bounds for LHS
		double yearsLimART = 2010;
		double yearsLimgen[2] = {2007, 2017 };
		double yearsLimmsm[2] = {1986, 2017};
		double yearsLimplha = 2015;
		double yearsLimpreg[2] = { 2006, 2017 };
		double yearsLimfsw = 2017;

		double artlimdatalow = 500 ;
		double artlimdataup = 2000 ;
		double prevmsmlimdatalow[2] = { 0.08, 0.08 };
		double prevmsmlimdataup[2] = { 0.6, 0.2 };
		double prevgenlimdatalow[2] = { 0.01, 0.008 };
		double prevgenlimdataup[2] = { 0.04, 0.025 };
		double prevpreglimdatalow[2] = { 0.004,0.005 };
		double prevpreglimdataup[2] = { 0.015, 0.013 };

		double plhalimdatalow = 3000;
		double plhalimdataup = 6000;
		double prevfswlimdatalow = 0.001;
		double prevfswlimdataup = 0.1;


#endif /* dataheader_h */
