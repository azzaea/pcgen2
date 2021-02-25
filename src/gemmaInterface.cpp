#include "gemma.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <Rcpp.h>
#define STRICT_R_HEADERS

using namespace Rcpp;
using namespace std;


// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar).

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/

// [[Rcpp::export]]
int rgemmamain(std::vector<std::string> strings)
{
  // int main(int argc, char *argv[])

  GEMMA cGemma;
  PARAM cPar;

  //std::vector<std::string> strings {"-g ", "/home/p287664/tmp/GEMMA-0.98.4/example/mouse_hs1940.geno.txt.gz",  "-p", "/home/p287664/tmp/GEMMA-0.98.4/example/mouse_hs1940.pheno.txt", "-a", "/home/p287664/tmp/GEMMA-0.98.4/example/mouse_hs1940.anno.txt", "-gk", "-o", "mouse_hs1940", "-nind",  "3"};

/* this works, but an alternative is also provided:
   std::vector<char*> cstrings;
  cstrings.reserve(strings.size());
  for(int i = 0; i < strings.size(); ++i){
    cstrings.push_back(const_cast<char*>(strings[i].c_str()));
  }
  */
std::vector<char*> cstrings{};

for(auto& string : strings)
{
  cstrings.push_back(&string.front());
}

  gsl_set_error_handler (&gemma_gsl_error_handler);

  if (cstrings.size() <= 1) {
    cGemma.PrintHeader();
    cGemma.PrintHelp(0);
    return EXIT_SUCCESS;
  }
/*  if (cstrings.size() == 2 && &cstrings[1][0] == '-' && &cstrings[1][1] == 'h') {
    cGemma.PrintHeader();
    cGemma.PrintHelp(0);
    return EXIT_SUCCESS;
  }
 if (cstrings.size() == 3 && inptargs[1][0] == '-' && inptargs[1][1] == 'h') {
    string str;
    str.assign(inptargs[2]);
    cGemma.PrintHeader();
    cGemma.PrintHelp(atoi(str.c_str()));
    return EXIT_SUCCESS;
  }
  if (cstrings.size() == 2 && inptargs[1][0] == '-' && inptargs[1][1] == 'l') {
    cGemma.PrintHeader();
    cGemma.PrintLicense();
    return EXIT_SUCCESS;
  }
*/
  //cGemma.PrintHeader();
  //cGemma.PrintLicense();
  //if(!cstrings.empty())
    cGemma.Assign(cstrings.size(), cstrings.data(), cPar);


  ifstream check_dir((cPar.path_out).c_str());
  if (!check_dir) {
#ifdef WINDOWS
    mkdir((cPar.path_out).c_str());
#else
    mkdir((cPar.path_out).c_str(), S_IRWXU | S_IRGRP | S_IROTH);
#endif
  }

  if (!is_quiet_mode())
    cGemma.PrintHeader();

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  if (is_quiet_mode()) {
    stringstream ss;
    cout.rdbuf(ss.rdbuf());
  }

  cPar.CheckParam();

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  cGemma.BatchRun(cPar);

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

 cGemma.WriteLog(cstrings.size(), cstrings.data(), cPar);

  return EXIT_SUCCESS;

//return true;
}


