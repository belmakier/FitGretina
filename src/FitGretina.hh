#ifndef FIT_GRETINA_HH
#define FIT_GRETINA_HH

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLatex.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include "Minuit2/Minuit2Minimizer.h"

class fit {
public:
  std::vector<std::string> names;
  std::vector<std::string> fnames;
  std::vector<std::string> histNames;
  std::vector<std::string> nondopHistNames;
  std::vector<TH1D*> hists;
  std::vector<double> scales;
  std::vector<double> low_errors;
  std::vector<double> up_errors;
  std::vector<std::string> types;

  std::vector<TH1D*> nondop_hists;

  TH1D *expthist;
  TH1D *nondop_expthist;
  TH1D *np_expthist;
  double np_scale;
  
  int rebin = 1;
  int evalct = 0;
  int binwidth = 1; //this represents the number of DTA bins going into the fitted histograms
  int method = 1; //1 = LL, 0 = chi2
  bool doubleExp;
  bool singleExp;
  bool scaledBkg;
  double nondop_chi2_scale;

  //background parameters
  double amp1, amp2;
  double tau1, tau2;
  TF1 *bkgd;

  //nondop background parameters
  double nondop_amp1, nondop_amp2;
  double nondop_tau1, nondop_tau2;
  TF1 *nondop_bkgd;

  double bkg_scale_slope;
  double bkg_scale_slope_nd;

  //fit range
  double low, high;

  fit(double l, double h);
  void load();  
  std::pair<double,double> eval(double x);
  std::pair<double,double> nondop_eval(double x);
  double operator()(const double *par);
  TGraph *getbgfunc(const double *par);
  TGraph *getfunc(const double *par);
  TGraph *nondop_getfunc(const double *par);
  TGraph *getfunc(int ind, const double *par);
  TGraph *nondop_getfunc(int ind, const double *par);
  TH1D* get_bs_hist();
  TH1D* get_nd_hist();
};

class fitManager {
public:
  fit *ff;

  double low;
  double high;
  std::string inpHistFileName;
  std::string histType; //1D/2D
  std::string dop_p_name;
  std::string dop_np_name;
  std::string nodop_name;
  int projx = -1;
  int gate = -1;
  double gate_low, gate_high;
  double np_scale = 1;
  std::string bkgType;
  bool doubleExp = true;
  bool singleExp = true;
  bool scaledBkg = false;
  bool nofit = false;
  double nondop_chi2_scale;
  std::string errorType;
  std::vector<std::string> ffiles;
  std::vector<std::string> labels;
  std::vector<std::string> histNames;
  std::vector<std::string> nondopHistNames;
  std::vector<std::string> types;
  std::vector<double> pars;

  std::vector<double> dopBgPars;
  std::vector<double> nodopBgPars;

  std::string outrootfile;
  std::string outdatfile;
  std::string outpdf;
  std::string ampfilename;

  double ROI_ll, ROI_ul;

  int rebin;
  std::string method;

  ROOT::Minuit2::Minuit2Minimizer *min=NULL;

  fitManager();

  void Read(std::string inpfile);
  void Setup();
  void PrintHeader(std::vector<std::string> names);
  void Load(std::string inpfile);
  double DoFit();
};

#endif
