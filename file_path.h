#pragma once
#ifndef FILES_PATH_h
#define FILES_PATH_h
#include <iostream>
#include <new>
#include <string>
using namespace std;

string *sigmafile, *bestrunfile, *postsizefile, *paramsfile, *runsoutfile, *runsoutfile_NA, *logpriorfile, *lhsfile, *out_file, *mean_file, *covmat_file;

void go_tofile() {

		postsizefile = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\postsize.txt" };
		paramsfile = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\mcmcsamples.txt" };
		runsoutfile = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\Runs.txt" };
		runsoutfile_NA = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\Runs_NA.txt" };
		sigmafile = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\paramsSigmas.txt" };
		bestrunfile = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\BestRuns.txt" };
		logpriorfile = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\logprior.txt" };
		lhsfile = new string{ "\\\\fi--didef2\\Tmp\\Juan\\BrazilWB\\LHSs\\LHS.txt" };
		out_file = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\MCMCchain.txt" };
		mean_file = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\means.txt" };
		covmat_file = new string{ "\\\\fi--didef2\\Tmp\\Juan\\Belize\\covmat.txt" };



}


#endif /* files_path_h */#pragma once
