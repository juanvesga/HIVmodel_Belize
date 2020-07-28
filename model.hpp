//
//  modelTBsear_hpp
//  TBsear
//
//  Created by juan fernando vesga on 01/03/2017.
//  Copyright © 2017 juan fernando vesga. All rights reserved.
//

#ifndef model_hpp
#define model_hpp
#include <stdio.h>
#include <vector>


//---------------------------------------------Class declaration--------------------------------------------------------------------------------
class Results {
    
public:
    Results(int Length); //constructor
    
    std::vector< std::vector<double>>state;
    
    ~Results(); //desctructor
    
    
private:
    int arraySize; //this is where we will store the size of the array.
    void CreateArray(); //custom function to create the array


};


class Model
{
public:
    Model();
    ~Model();
    //
    Results *  Run(double * ImportedParams, int sizepars , int t0, int tend , std::vector<double> interv);
};


//---------------------------------------------Constant Model parameters --------------------------------------------------------------------------------
// Use defintions to index through X states
// Define Stage of disease (S) 
#define Sus 0
#define Ac 1
#define U_500p 2
#define U_500 3
#define U_350 4
#define U_200 5
#define D_500p 6
#define D_500 7
#define D_350 8
#define D_200 9
#define A_500p 10
#define A_500 11
#define A_350 12
#define A_200 13
#define S_500p 14
#define S_500 15
#define S_350 16
#define S_200 17
#define Di_500p 18
#define Di_500 19
#define Di_350 20
#define Di_200 21

// Define Risk Group (R)
#define nr_w 0  //no risk
#define lo_w 1  
#define fsw 2
#define nr_m 3
#define lo_m 4
#define cfsw 5
#define bi 6
#define msm 7


//Define outputs
#define prev 704                                                                         //---------->
#define prevlo_w 705
#define prevfsw 706
#define prevlo_m 707
#define prevcfsw 708
#define prevbi 709
#define prevmsm 710
#define prevpreg15 711//+ 9
#define prevpreg34 712//+ 10
#define prevpreg49 713//+ 11
#define prevM15 714//+ 12
#define prevM34 715//+ 13
#define inc  716//+ 14                                                                          //---------->
#define inclo_w 717//+ 15
#define incfsw 718//+ 16
#define inclo_m 719//+ 18
#define inccfsw 720//+ 19
#define incbi 721//+ 20
#define incmsm 722//+ 21
#define sus  723//+ 23                                                                           //---------->
#define suslo_w 724//+ 24
#define susfsw 725//+ 25
#define suslo_m 726//+ 27
#define suscfsw 727//+ 28
#define susbi 728//+ 29
#define susmsm 729//+ 30
#define undet 730//+ 32                                                                          //---------->
#define detec 731//+ 33                                                                          //---------->
#define onart 732//+ 34                                                                          //---------->
#define suppr 733//+ 35                                                                          //---------->
#define diseng 734//+ 36                                                                          //---------->
#define mortundet 735//+ 37                                                                          //---------->
#define mortdetec 736//+ 38                                                                          //---------->
#define mortonart 737//+ 39                                                                          //---------->
#define mortsuppr 738//+ 40                                                                          //---------->
#define mortdiseng 739//+ 41                                                                          //---------->
#define mortaHIV 740//+ 42                                                                          //---------->
#define mortaAIDS 741//+ 43                                                                          //---------->
#define newdetec 742//+ 44                                                                          //---------->
#define newart 743//+ 45                                                                          //---------->
#define dxdeath 744//+ 46                                                                          //---------->
#define dxaids 745//+ 47                                                                          //---------->
#define dxhiv 746//+ 48                                                                          //---------->
#define lifeyears 747//+ 49                                                                          //---------->
#define plha 748//+ 50                                                                          //---------->
#define pop 749//+ 51                                                                          //---------->
#define poplo_w 750//+ 52
#define popfsw 751//+ 53
#define poplo_m 752//+ 55
#define popcfsw 753//+ 56
#define popbi 754//+ 57
#define popmsm 755//+ 58
#define poppreg15 756//+ 60
#define poppreg34 757//+ 61
#define poppreg49 758//+ 62
#define popM15 759//+ 63
#define popM34 760//+ 64
#define newinf 761                                                                          //---------->
#define newinflo_w 762//+ 66
#define newinffsw 763//+ 67
#define newinflo_m 764//+ 69
#define newinfcfsw 765//+ 70
#define newinfbi 766//+ 71
#define newinfmsm 767//+ 72
#define prevmsmall 768//+ 77
#define incmsmall 769//+ 78
#define newinfmsmall 770//+ 79
#define notaids_a0 771//+ 80
#define notaids_a1 772//+ 81
#define notaids_a2 773//+ 82
#define notaids_a3 774//+ 83
#define nothiv_a0 775//+ 84
#define nothiv_a1 776//+ 85
#define nothiv_a2 777//+ 86
#define nothiv_a3 778//+ 87
#define notdead_a0 779//+ 88
#define notdead_a1 780//+ 89
#define notdead_a2 781//+ 90
#define notdead_a3 782//+ 91


/////<<<<<<<<<<<ALL MODELS PARAMETERS<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

//-------->>>>>>>>>>>>>>	 =
//	 =	SIMULATION	CONDITIONS
//-------->>>>>>>>>>>>>>>	 =
const int nstages = (22 * 8 * 4);
const int S_all = (22*8*4) + 79; // 
const int S = 22; // Model stages
const int R = 8; //Risk groups
const int A = 4; // Age-years
const double dt = 0.1;// //time step
const int timeZero = 1975;
const int noutcomes = 79;
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//:::::::::::::::::::::::::::::: FIXED PARAMETERS - for all regions - 
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//Test
const double test_start=1984;
const double hivtest_start = 2000;
const double AIDStestRR = 2.5;
const double  Undetected_Fraction = 0.8;//UndFrac
const double undeR = 1; 
// Population Structure
const double f_prop = 0.5; 
const double m_prop = 1-0.5;
const double NoRM_frac = 0.15;//	
const double NoRW_frac = 0.2;
const double Bi_frac = 0.5;
const double turn_SW = 0.06666;//

const double ageDist[] = { 0.352115988,0.418218595,0.148499175,0.081166242};
const double initPop =  70140; // POP 15 to 100 in belize 1975

const double Tratios[R][A] = {
{   0.701996928,	1,           	0.655913978,	0.322580645},
{   0.701996928,	1,	            0.655913978,	0.322580645 },
{	0.701996928,	1,           	0.655913978,	0.322580645 },
{	0.153894246,	0.559139785,	0.675399542,	0.656410256 },
{	0.153894246,	0.559139785,	0.675399542,	0.656410256 },
{	0.153894246,	0.559139785,	0.675399542,	0.656410256 },
{	0.153894246,	0.559139785,	0.675399542,	0.656410256 },
{	0.153894246,	0.559139785,	0.675399542,	0.656410256 }
};


// HIV Natural History 
const double durAc = 0.3;
//const double Acute_coef = 9.17;
const double lateStageHR = 0.275;// Hollingsworth
const double dur500p[2] = { 1/4.87 ,1/5.26};// Mangal (15.91)    male    Female
const double dur500[2] =  { 1/2.94 ,1/3.17};
const double dur350[2] =  { 1/4.76 ,1/5.14};
const double dur200[2] =  { 1/2.79 ,1/3.02};
const double mortHIV[4] = { 0.004110107,	0.011670254,	0.012889907,	0.344};// According to mangal etal 
const double surv_ARText = 4;// reduction calculated as 0.344/0.086 which is mortality at 200 as above over mortlity after 1yr ARt as estimated by Tuboi 2009 JAIDS
const double disRR = 1.5;//increased mort by disengagment assump
const double postAcute[4] = { 0.5251,0.2315,0.2401,0.0033 };// PostAcute progression (proportion 200-350-500-500p) Mangal et al
// Condoms

const double cond_msm = 0.75;
const double cond_fsw=0.8;
const double cond_het =  0.45;
const double cond_start=1984;
const double Condom_eff = 0.95; // condom effectiveness

// Mixing 
const double theta[R]={0.5 , 0.5, 0.8 ,0.5,0.5,0.5,0.5,0.5};
//const double propBiacts = 0.1; 
//const double c_f = 1.5; 
//const double c_m = 2.5;
//const double c_msm = 10;;// 34;
//const double c_bi = 11;// 35;
//const double c_cfsw = 8;// 27;
//const double c_fsw = 170;
//const double	ratio15_25 = 0.65;// 1.0000;	//Reference:ratio	of	casual partners	pryr	by	age
//const double	ratio25_50 = 1;// 0.625;	//PCAP	2008
//const double	ratio50_65 =  0.473;	//PCAP	2008
//const double	ratio65_99 =  0.36;	//PCAP	2008

const double	ratio15_25 = 1;// 3;// 1.0000;	//Reference:ratio	of	casual partners	pryr	by	age
const double	ratio25_50 = 0.75;// 1.5;// 0.625;	//PCAP	2008
const double	ratio50_65 = 0.50;	//PCAP	2008
const double	ratio65_99 = 0.3;// 0.2;	//PCAP	2008


const double    actsSW = 2;
const double	actsHet = 80.0;
//const double	actsMSM = 10;
const double	actsBi = 10;
const double	actsCFSW = 2;
// Interventions
const double supp = 0.85; // rate from ARt to suppressed
const double dis = 0.11; // Drop-ut rate
const double backRR = 1.1; // RR for fast rate from drop-out to ART




#endif /* model_new_hpp */
