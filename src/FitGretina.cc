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

#include "FitGretina.hh"

fit::fit(double l, double h) {
  low = l;
  high = h;
  doubleExp = true;
  singleExp = true;
  bkgd = new TF1("bkgd", "[0]*exp(-x/[1]) + [2]*exp(-x/[3])", low, high);
  nondop_bkgd = new TF1("nondop_bkgd", "[0]*exp(-x/[1]) + [2]*exp(-x/[3])", low, high);
  nondop_chi2_scale = 1.0;
}

void fit::load() {
  for (int i=0; i<fnames.size(); ++i) {
    TFile *file = new TFile(fnames[i].c_str());
    hists.push_back((TH1D*)file->Get<TH1F>(histNames[i].c_str()));
        
    hists[i]->SetDirectory(0);
    hists[i]->Rebin(rebin);
    nondop_hists.push_back((TH1D*)file->Get<TH1F>(nondopHistNames[i].c_str()));
    nondop_hists[i]->SetDirectory(0);
    nondop_hists[i]->Rebin(rebin);
    file->Close();
    delete file;
  }
}
  
std::pair<double,double> fit::eval(double x) {
  std::pair<double,double> val;
  if (!scaledBkg) {
    val.first = bkgd->Eval(x);
  }
  else {
    val.first = 0;
  }
  val.second = 0;
  for (int i=0; i<hists.size(); ++i) {
    double cts = hists[i]->GetBinContent(hists[i]->GetXaxis()->FindBin(x));
    if (!types[i].compare("BKG")) {
      cts = cts * (1.0 + (x-500.0)*bkg_scale_slope);
    }
    val.first += scales[i]*cts*binwidth;
    val.second += cts*scales[i]*scales[i]*binwidth*binwidth;
  }
  val.second = std::sqrt(val.second);
  return val;
}

std::pair<double,double> fit::nondop_eval(double x) {
  std::pair<double,double> val;
  //double val = nondop_bkgd->Eval(x);
  if (!scaledBkg) {
    val.first = nondop_bkgd->Eval(x);
  }
  else {
    val.first = 0.0;
  }
  val.second = 0;
  for (int i=0; i<nondop_hists.size(); ++i) {
    //val += scales[i]*nondop_hists[i]->GetBinContent(nondop_hists[i]->GetXaxis()->FindBin(x));
    double cts = nondop_hists[i]->GetBinContent(nondop_hists[i]->GetXaxis()->FindBin(x));
    if (!types[i].compare("BKG")) {
      cts = cts * (1.0 + (x-500.0)*bkg_scale_slope_nd);
    }
    val.first += scales[i]*cts*binwidth;
    val.second += cts*scales[i]*scales[i]*binwidth*binwidth;
  }
  val.second = std::sqrt(val.second);
  return val;
}
  
double fit::operator()(const double *par) {
  if (!scaledBkg) {
    amp1 = par[0];
    tau1 = par[1];
    amp2 = par[2];
    tau2 = par[3];

    nondop_amp1 = par[4];
    nondop_tau1 = par[5];
    nondop_amp2 = par[6];
    nondop_tau2 = par[7];

    for (int i=0; i<4; ++i) {
      bkgd->SetParameter(i, par[i]);
    }

    for (int i=4; i<8; ++i) {
      nondop_bkgd->SetParameter(i-4, par[i]);
    }
  }
  else {
    bkg_scale_slope = par[0];
    bkg_scale_slope_nd = par[1];
  }

  for (int i=0; i<scales.size(); ++i) {
    scales[i] = par[8+i];
  }

  double chi2 = 0;
  if (method==0) {
    for (int i=1; i<=expthist->GetNbinsX(); ++i) {
      if (expthist->GetBinCenter(i) < low || expthist->GetBinCenter(i) > high) { continue; }
      std::pair<double, double> val;
      val = eval(expthist->GetBinCenter(i));
      if (expthist->GetBinContent(i) > 0) {
        chi2 += std::pow((expthist->GetBinContent(i) - val.first), 2)/(std::pow(expthist->GetBinError(i), 2) + std::pow(val.second, 2));
      }
      else {
        chi2 += std::pow((expthist->GetBinContent(i) - val.first), 2)/(1.0 + std::pow(val.second, 2));
      }
    }
    for (int i=1; i<=nondop_expthist->GetNbinsX(); ++i) {  
      if (nondop_expthist->GetBinCenter(i) < low || nondop_expthist->GetBinCenter(i) > high) { continue; }
      std::pair<double, double> val;
      val = nondop_eval(expthist->GetBinCenter(i));
      if (nondop_expthist->GetBinContent(i) > 0) {
        chi2 += nondop_chi2_scale * std::pow((nondop_expthist->GetBinContent(i) - val.first), 2)/(std::pow(nondop_expthist->GetBinError(i), 2) + std::pow(val.second, 2));
      }
      else {
        chi2 += nondop_chi2_scale * std::pow((nondop_expthist->GetBinContent(i) - val.first), 2)/(1.0 + std::pow(val.second, 2));
      }
    }
  }
  else if (method == 1) {
    //log likelihood
    double nll = 0;
    for (int i=1; i<=expthist->GetNbinsX(); ++i) {
      if (expthist->GetBinCenter(i) < low || expthist->GetBinCenter(i) > high) { continue; }
      std::pair<double, double> val = eval(expthist->GetBinCenter(i));
      val.first += np_scale * np_expthist->GetBinContent(i);
      double y = expthist->GetBinContent(i);
      if (y > 0) {
        val.first = std::max(val.first, 1e-5);
        nll += (val.first + y * TMath::Log(y/val.first) - y);
      }
      else {
        nll += val.first;
      }
    }

    for (int i=1; i<=nondop_expthist->GetNbinsX(); ++i) {
      if (nondop_expthist->GetBinCenter(i) < low || nondop_expthist->GetBinCenter(i) > high) { continue; }
      std::pair<double, double> val = nondop_eval(nondop_expthist->GetBinCenter(i));
      double y = nondop_expthist->GetBinContent(i);
      if (y > 0) {
        val.first = std::max(val.first, 1e-5);
        nll += nondop_chi2_scale * (val.first + y * TMath::Log(y/val.first) - y);
      }
      else {
        nll += nondop_chi2_scale * val.first;
      }
    }
    chi2 = 2.0*nll;
  }

  if (evalct %100 == 0) {
    if (doubleExp) {
      std::cout << "\r" << Form("%9.8g", chi2) << " : ";
      for (int i=0; i<std::min((int)(8+scales.size()), 15); ++i) {
        std::cout << Form("%9.3g  ", par[i]);
      }
      std::cout << std::flush;
    }
    else if (singleExp) {
      std::cout << "\r" << Form("%9.8g", chi2) << " : ";
      for (int i=0; i<std::min((int)(8+scales.size()), 17); ++i) {
        if (i==2 || i==3 || i==6 || i==7) { continue; }
        std::cout << Form("%9.3g  ", par[i]);
      }
      std::cout << std::flush;
    }
    else {
      std::cout << "\r" << Form("%9.8g", chi2) << " : ";
      for (int i=0; i<std::min((int)(8+scales.size()), 17); ++i) {
        if (i<=7) { continue; }
        std::cout << Form("%9.3g  ", par[i]);
      }
      std::cout << std::flush;

    }
  }
  evalct += 1;
  return chi2;    
}

TGraph *fit::getbgfunc(const double *par) {
  double chi2 = (*this)(par); //this updates the internal values
  TGraph *gr = new TGraph();
  gr->SetLineColor(kGreen+2);
  for (int i=1; i<=expthist->GetNbinsX(); ++i) {
    double val = 0;
    for (int ind=0; ind <scales.size(); ++ind) {
      if (!types[ind].compare("BKG")) {
        val += scales[ind]*hists[ind]->GetBinContent(i)*binwidth;
      }
    }
    gr->AddPoint(expthist->GetXaxis()->GetBinCenter(i), val);
  }
  return gr;
}

TGraph *fit::getfunc(const double *par) {
  double chi2 = (*this)(par); //this updates the internal values
  TGraph *gr = new TGraph();
  gr->SetLineColor(kRed);
  for (int i=1; i<=expthist->GetNbinsX(); ++i) {
    gr->AddPoint(expthist->GetBinCenter(i), eval(expthist->GetBinCenter(i)).first);
  }
  return gr;
}

TGraph *fit::nondop_getfunc(const double *par) {
  double chi2 = (*this)(par); //this updates the internal values
  TGraph *gr = new TGraph();
  gr->SetLineColor(kRed);
  for (int i=1; i<=nondop_expthist->GetNbinsX(); ++i) {
    gr->AddPoint(nondop_expthist->GetBinCenter(i), nondop_eval(nondop_expthist->GetBinCenter(i)).first);
  }
  return gr;
}

  
TGraph *fit::getfunc(int ind, const double *par) {
  double chi2 = (*this)(par); //this updates the internal values
  TGraph *gr = new TGraph();
  gr->SetLineColor(kRed);
  for (int i=1; i<=expthist->GetNbinsX(); ++i) {
    gr->AddPoint(expthist->GetBinCenter(i), binwidth*scales[ind]*hists[ind]->GetBinContent(hists[ind]->GetXaxis()->FindBin(expthist->GetBinCenter(i))));
  }
  return gr;
}

TGraph *fit::nondop_getfunc(int ind, const double *par) {
  double chi2 = (*this)(par); //this updates the internal values
  TGraph *gr = new TGraph();
  gr->SetLineColor(kRed);
  for (int i=1; i<=nondop_expthist->GetNbinsX(); ++i) {
    gr->AddPoint(nondop_expthist->GetBinCenter(i), binwidth*scales[ind]*nondop_hists[ind]->GetBinContent(nondop_hists[ind]->GetXaxis()->FindBin(nondop_expthist->GetBinCenter(i))));
  }
  return gr;
}

TH1D* fit::get_bs_hist() {
  TH1D* expthist_bs = (TH1D*)expthist->Clone("expthist_bs");
  expthist_bs->Add(np_expthist, -np_scale);
  return expthist_bs;
}

TH1D* fit::get_nd_hist() {
  TH1D* expthist_nd = (TH1D*)nondop_expthist->Clone("expthist_nd");
  return expthist_nd;
}

fitManager::fitManager() {
  method="LOGLIKELIHOOD";
  histType = "2D";
}

void fitManager::Load(std::string inpfile) {
  Read(inpfile);
  Setup();
}  

void fitManager::Read(std::string inpfile) {
  FILE *file;
  if ( !(file = fopen(inpfile.c_str(), "ra") ) ) {
    std::cerr << "Warning! Input file " << inpfile << " not found!" << std::endl;      
    exit(1);
  }

  std::stringstream ss;
  char cline[2048];

  rebin = 1;

  while(std::fgets(cline, sizeof cline, file)!=NULL) {    
    std::string line(cline);
    if (line.size() == 0) { continue; }
    if (line[0] == '#') { continue; }
    if (line[0] == ';') { continue; }

    ss.clear();
    ss.str(line);

    std::string keyword;
    ss >> keyword;

    if (!keyword.compare("INPUT")) {
      ss >> inpHistFileName;
    }
    else if (!keyword.compare("ID")) {
      std::cout << "ID keyword depreciated!" << std::endl;
      //ss >> id;
    }
    else if (!keyword.compare("COMB")) {
      std::cout << "COMB keyword depreciated!" << std::endl;
      //ss >> comb;
    }
    else if (!keyword.compare("DTAGATE")) {
      std::cout << "DTAGATE keyword depreciated!" << std::endl;
      //ss >> dta_low >> dta_high;
    }
    else if (!keyword.compare("PPARGATE")) {
      std::cout << "DTAGATE keyword depreciated!" << std::endl;
      //ss >> ppar_low >> ppar_high;
    }
    else if (!keyword.compare("TYPE")) {
      ss >> histType;
    }
    else if (!keyword.compare("DOP_P")) {
      ss >> dop_p_name;
    }
    else if (!keyword.compare("DOP_NP")) {
      ss >> dop_np_name;
    }
    else if (!keyword.compare("NODOP")) {
      ss >> nodop_name;
    }
    else if (!keyword.compare("PROJX")) {
      ss >> projx;   
    }
    else if (!keyword.compare("GATE")) {
      gate=1;
      ss >> gate_low >> gate_high;
    }
    else if (!keyword.compare("NP_SCALE")) {
      ss >> np_scale;
    }
    else if (!keyword.compare("ERROR")) {
      ss >> errorType;
    }
    else if (!keyword.compare("ROI")) {
      ss >> ROI_ll >> ROI_ul;
    }
    else if (!keyword.compare("FFILE")) {
      std::string type;
      std::string fname;
      std::string hname;
      std::string nd_hname;
      std::string label;
      double par;
      ss >> type >> fname >> hname >> nd_hname >> label >> par;
      types.push_back(type);
      ffiles.push_back(fname);
      histNames.push_back(hname);
      nondopHistNames.push_back(nd_hname);
      labels.push_back(label);
      pars.push_back(par);
    }
    else if (!keyword.compare("DOPBKG")) {
      double a,b,c,d;
      ss >> a >> b >> c >> d;
      dopBgPars.push_back(a);
      dopBgPars.push_back(b);
      dopBgPars.push_back(c);
      dopBgPars.push_back(d);
    }
    else if (!keyword.compare("NDBKG")) {
      double a,b,c,d;
      ss >> a >> b >> c >> d;
      nodopBgPars.push_back(a);
      nodopBgPars.push_back(b);
      nodopBgPars.push_back(c);
      nodopBgPars.push_back(d);
    }
    else if (!keyword.compare("RANGE")) {
      ss >> low >> high;
    }
    else if (!keyword.compare("OUTDAT")) {
      ss >> outdatfile;
    }
    else if (!keyword.compare("OUTROOT")) {
      ss >> outrootfile;
    }
    else if (!keyword.compare("AMPFILE")) {
      ss >> ampfilename;
    }
    else if (!keyword.compare("OUTPDF")) {
      ss >> outpdf;
    }
    else if (!keyword.compare("REBIN")) {
      ss >> rebin;
    }
    else if (!keyword.compare("METHOD")) {
      ss >> method;
    }
    else if (!keyword.compare("BKGTYPE")) {
      ss >> bkgType;
      if (!bkgType.compare("EXP")) {
        singleExp = true;
        doubleExp = false;
        scaledBkg = false;
      }
      else if (!bkgType.compare("DOUBLEEXP")) {
        singleExp = true;
        doubleExp = true; 
        scaledBkg = false;
      }
      else if (!bkgType.compare("SCALED")) {
        singleExp = false;
        doubleExp = false; 
        scaledBkg = true;
      }
      else if (!bkgType.compare("NONE")) {
        singleExp = false;
        doubleExp = false; 
        scaledBkg = false;
      }
    }
    else if (!keyword.compare("ND_CHI2_SCALE")) {
      ss >> nondop_chi2_scale;
    }
    else if (!keyword.compare("NOFIT")) {
      nofit = true;
    }
  }
  //check if options are valid
  if (dop_p_name.size() == 0 || dop_np_name.size()==0 || nodop_name.size()==0){
    std::cerr << "Prompt, non-prompt, and non-Doppler corrected histograms must be specified with DOP_P, DOP_NP, and NODOP!" << std::endl;
  }
  if (histType.compare("1D") && histType.compare("2D")) {
    std::cerr << "Valid TYPE options are 1D and 2D!" << std::endl;
  }
  if (ROI_ll == ROI_ul) { 
    std::cerr << "Specify region of interest for chi^2 printout with ROI!" << std::endl;
  }
  if (gate_low==gate_high && projx==-1 && !histType.compare("2D")) {
    std::cerr << "Either GATE or PROJX option must be specified!" << std::endl;
  }
  if (gate>0 && projx>0) {
    std::cerr << "Specify only one of GATE or PROJX!" << std::endl;
  }
  if (low==high) {
    std::cerr << "Specify fit range with RANGE!" << std::endl;
  }
}

void fitManager::Setup() {
  ff = new fit(low,high);
  ff->rebin = rebin;
  if (!method.compare("CHI2")) {
    ff->method = 0;
  }
  else if (!method.compare("LOGLIKELIHOOD")) {
    ff->method = 1;
  }

  TFile *expt = new TFile(inpHistFileName.c_str());
  ff->np_scale = np_scale;
  ff->nondop_chi2_scale = nondop_chi2_scale;
  if (!histType.compare("2D")) {
    TH2F *expt_h2 = expt->Get<TH2F>(dop_p_name.c_str());
    TH2F *expt_h2_np = expt->Get<TH2F>(dop_np_name.c_str());
    TH2F *nondop_expt_h2 = expt->Get<TH2F>(nodop_name.c_str());

    expt_h2->SetDirectory(0);
    expt_h2_np->SetDirectory(0);
    nondop_expt_h2->SetDirectory(0);

    if (projx>0) {
      ff->expthist = expt_h2->ProjectionX("expthist", projx, projx);
      ff->np_expthist = expt_h2_np->ProjectionX("expthist_np", projx, projx);
      ff->nondop_expthist = nondop_expt_h2->ProjectionX("nondop_expthist", projx, projx);
    }
    else if (gate) {
      int binlo = expt_h2->GetYaxis()->FindBin(gate_low);
      int binhi = expt_h2->GetYaxis()->FindBin(gate_high) - 1;
      std::cout << "Gate = (" << gate_low << ", " << gate_high << ")" << std::endl;

      ff->binwidth = binhi-binlo + 1;

      ff->expthist = expt_h2->ProjectionX("expthist", binlo, binhi); 
      ff->np_expthist = expt_h2_np->ProjectionX("expthist_np", binlo, binhi);
      ff->nondop_expthist = nondop_expt_h2->ProjectionX("nondop_expthist", binlo, binhi);      
    }
    else { 
      std::cerr << "Invalid gate/projection options for 2D histograms! Incomplete setup!" << std::endl;
      exit(1);
    }
  }
  else if (!histType.compare("1D")) {
    ff->expthist = (TH1D*)expt->Get<TH1F>(dop_p_name.c_str());
    ff->np_expthist = (TH1D*)expt->Get<TH1F>(dop_np_name.c_str());
    ff->nondop_expthist = (TH1D*)expt->Get<TH1F>(nodop_name.c_str());
  }
  else {
    std::cerr << "Invalid hist type option!" << std::endl;
    exit(1);
  }

  
  ff->expthist->SetDirectory(0);
  ff->np_expthist->SetDirectory(0);
  ff->nondop_expthist->SetDirectory(0);

  ff->expthist->RebinX(rebin);
  ff->np_expthist->RebinX(rebin);
  ff->nondop_expthist->RebinX(rebin);

  for (int i=0; i<ffiles.size(); ++i) {
    ff->scales.push_back(pars[i]);
    ff->low_errors.push_back(0.0);
    ff->up_errors.push_back(0.0);
    ff->fnames.push_back(ffiles[i]);
    ff->names.push_back(labels[i]);
    ff->histNames.push_back(histNames[i]);
    ff->nondopHistNames.push_back(nondopHistNames[i]);
    ff->types.push_back(types[i]);
  }

  expt->Close();
  delete expt;
  ff->load();
}

void fitManager::PrintHeader(std::vector<std::string> names) {
  if (doubleExp) {
    std::cout << "chi2     " << " : ";
    for (int i=0; i<std::min((int)(names.size()), 15); ++i) {
      std::cout << Form("%9s  ", names[i].substr(0, std::min(9, (int)names[i].size())).c_str());
    }
    std::cout << std::endl;
  }
  else if (singleExp){
    std::cout << "chi2     " << " : ";
    for (int i=0; i<std::min((int)(names.size()), 17); ++i) {
      if (i==2 || i==3 || i==6 || i==7) { continue; }
      std::cout << Form("%9s  ", names[i].substr(0, std::min(9, (int)names[i].size())).c_str());
    }
    std::cout << std::endl;
  }
  else if (scaledBkg) {
    std::cout << "chi2     " << " : ";
    for (int i=0; i<std::min((int)(names.size()), 17); ++i) {
      if (i>=2 && i<=7) { continue; }
      std::cout << Form("%9s  ", names[i].substr(0, std::min(9, (int)names[i].size())).c_str());
    }
    std::cout << std::endl;
  }
  else {
    std::cout << "chi2     " << " : ";
    for (int i=0; i<std::min((int)(names.size()), 17); ++i) {
      if (i<=7) { continue; }
      std::cout << Form("%9s  ", names[i].substr(0, std::min(9, (int)names[i].size())).c_str());
    }
    std::cout << std::endl;
  }
}

double fitManager::DoFit() {
  pars.clear();
  
  for (int i=0; i<4; ++i) {
    pars.push_back(dopBgPars[i]);
  }
  for (int i=0; i<4; ++i) {
    pars.push_back(nodopBgPars[i]);
  }
  std::vector<std::string> names = {"amp1", "tau1", "amp2", "tau2",
    "nd_amp1", "nd_tau1", "nd_amp2", "nd_tau2"};
  std::vector<double> LL = {0.0, 200.0, 0.0, 200.0, 0.0, 200.0, 0.0, 200.0};
  std::vector<double> UL = {1e4,5e3,1e4,5e3,1e4,5e3,1e4,5e3};
  if (scaledBkg) {
    names[0] = "bkg_slope";
    names[1] = "bkg_slope_nd";
    pars[0] = dopBgPars[0];
    pars[1] = nodopBgPars[0];
    LL[0] = -0.1;
    LL[1] = -0.1;
    UL[0] = 0.1;
    UL[1] = 0.1;
  }
  for (int i=0; i<ff->scales.size(); ++i) {
    pars.push_back(ff->scales[i]);
    names.push_back(ff->names[i]);
    LL.push_back(0.0);
    UL.push_back(1.0);
  }

  //create minimizer
  if (min != NULL) { delete min; }
  
  min = new ROOT::Minuit2::Minuit2Minimizer("migrad");

  //ROOT::Minuit2::Minuit2Minimizer *min=new ROOT::Minuit2::Minuit2Minimizer("simplex");
  
  if (min==NULL) { std::cerr << "Minimizer failed to create!" << std::endl; exit(1); }

  ff->doubleExp = doubleExp;
  ff->singleExp = singleExp;
  ff->scaledBkg = scaledBkg;

  ROOT::Math::Functor f_init(*ff,ff->scales.size() + 8);
	
  min->SetErrorDef(1.);
  min->SetMaxFunctionCalls(100000);
  min->SetMaxIterations(100000);
  //min->SetTolerance(0.001);
  //min->SetPrecision(1e-8);
  min->SetFunction(f_init);
	
  for( unsigned int i = 0; i < pars.size(); ++i ) {
    min->SetLimitedVariable(i, names.at(i), pars.at(i), 0.01, LL.at(i), UL.at(i));
    if (i <=7) {      
      min->SetVariableStepSize(i, 1.0);
    }
    else {
      min->SetVariableStepSize(i, pars.at(i)/100.0);
    }
  }

  if (scaledBkg) {
    min->SetFixedVariable(2, names[2], 0);
    min->SetFixedVariable(3, names[3], 1.);
    min->SetFixedVariable(4, names[3], 0.);
    min->SetFixedVariable(5, names[3], 1.);
    min->SetFixedVariable(6, names[6], 0);
    min->SetFixedVariable(7, names[7], 1.);
  }
  if (!doubleExp && !scaledBkg) {
    min->SetFixedVariable(2, names[2], 0);
    min->SetFixedVariable(3, names[3], 1.);
    min->SetFixedVariable(6, names[6], 0);
    min->SetFixedVariable(7, names[7], 1.);
  }
  if (!singleExp && !scaledBkg) {
    min->SetFixedVariable(0, names[0], 0);
    min->SetFixedVariable(1, names[1], 1.);
    min->SetFixedVariable(4, names[4], 0);
    min->SetFixedVariable(5, names[5], 1.);
  }

  PrintHeader(names);	
  min->FixVariable(1);
  min->FixVariable(5);
  if (nofit) { //fix everything
    for (int i=0; i<pars.size(); ++i) {
      min->FixVariable(i);
    }
  }
  min->Minimize();

  std::cout << std::endl;

  min->PrintResults();
  min->ReleaseVariable(1);
  min->ReleaseVariable(5);
  
  PrintHeader(names);	
  if (!nofit) {
    min->Minimize();
  }

  std::cout << std::endl;
  
  min->PrintResults();
  //min->SetMinimizerType(ROOT::Minuit2::EMinimizerType::kCombined);
  if (!nofit) {
    int ct = 5;
    while (ct > 0) {
      if (min->Status() != 0) {
        for( unsigned int i = 0; i < pars.size(); ++i ) {
          if (i<2 && scaledBkg) { continue; }
          double ul = std::max(0.0001, 3.0*min->X()[i]);
          min->SetLimitedVariable(i, names.at(i), min->X()[i], 0.01, LL.at(i), ul);
          UL[i] = ul;
          min->SetVariableStepSize(i, min->X()[i]/100.0);
        }
        PrintHeader(names);	
        min->Minimize();
        std::cout << std::endl;
        min->PrintResults();
        std::cout << std::endl;
        ct -= 1;
      }
      else { break; }
    }
  }

  /*
    std::cout << "MINOS ERRORS" << std::endl;
    for (int i=0; i<names.size(); ++i) {
    double low, high;
    min->GetMinosError(i, low, high);
    std::cout << i << "  " << low << "  " << high << std::endl;
    }
  */
	
  min->PrintResults();

  if (!errorType.compare("MINOS")) {
    std::cout << "MINOS ERRORS" << std::endl;
  }
  std::ofstream ampfile(ampfilename);
  PrintHeader(names);
  for (int i=0; i<names.size(); ++i) {
    if (!errorType.compare("MINOS") && !nofit) {
      double low, high;
      if (!singleExp && (i<=7)) { 
        ampfile << i << "   " << names.at(i) << "   " << 0.0 << "   " << 0.0 << std::endl;
        continue;
      }
      if (!doubleExp && (i==2 || i==3 || i==6 || i==7)) { 
        ampfile << i << "   " << names.at(i) << "   " << 0.0 << "   " << 0.0 << std::endl;
        continue;
      }
      bool retval = min->GetMinosError(i, low, high);
      std::cout << std::endl;
      if (!retval) {
        std::cout << "doing my own errors" << std::endl;
        ROOT::Minuit2::Minuit2Minimizer *new_min = new ROOT::Minuit2::Minuit2Minimizer();

        new_min->SetErrorDef(1.);
        new_min->SetMaxFunctionCalls(100000);
        new_min->SetMaxIterations(100000);
        new_min->SetFunction(f_init);

        for( unsigned int i = 0; i < pars.size(); ++i ) {
          if (i <=7) {      
            new_min->SetLimitedVariable(i, names.at(i), min->X()[i], 1.0, LL.at(i), UL.at(i));
          }
          else {
            new_min->SetLimitedVariable(i, names.at(i), min->X()[i], 0.01, LL.at(i), UL.at(i));
            new_min->SetVariableStepSize(i, min->X()[i]/100.0);
          }
        }

        if (!doubleExp) {
          new_min->SetFixedVariable(2, names[2], 0);
          new_min->SetFixedVariable(3, names[3], 1.);
          new_min->SetFixedVariable(6, names[6], 0);
          new_min->SetFixedVariable(7, names[7], 1.);
        }
        if (!singleExp) {
          new_min->SetFixedVariable(0, names[0], 0);
          new_min->SetFixedVariable(1, names[1], 1.);
          new_min->SetFixedVariable(4, names[4], 0);
          new_min->SetFixedVariable(5, names[5], 1.);
        }

        new_min->Minimize();
        double opt_val = min->X()[i];
        double min_val = std::max(0.0, min->X()[i]-3.0*min->Errors()[i]);
        double max_val = std::max(0.00001, min->X()[i]+3.0*min->Errors()[i]);
        if (min->Errors()[i]<=1e-6) {
          min_val = 0.0;
          max_val = std::max(0.00001, opt_val*2.0);
        }
        max_val = std::min(max_val, UL.at(i));
        min_val = std::max(min_val, LL.at(i));
        std::cout << std::endl << "   " << opt_val << "  " << min_val << "  " << max_val << std::endl;
        double min_chi2 = min->MinValue();
        double val = (max_val + min->X()[i])/2.;
        double stepsize = (max_val - opt_val)/2.;
        int last_dir = 0;
        //upper error
        int ct = 100;
        int zeroct = 0;
        double last_gt_chi2 = -1;
        double last_lt_chi2 = min_chi2;
        double last_gt_val = -1;
        double last_lt_val = opt_val;
        std::cout << "finding upper limit" << std::endl;
        while (ct > 0) {
          new_min->SetFixedVariable(i, names[i], val);
          new_min->Minimize();
          //min->PrintResults();
          std::cout << std::endl << opt_val << "  " << val << "   " << min_val << "  " << max_val << "  " << new_min->MinValue() << "  " << min_chi2 << std::endl;
          
          if (std::abs(new_min->MinValue() - (min_chi2 + 1)) < 0.001) { break; }
          else if (new_min->MinValue() < min_chi2 - 1.0) { 
            std::cout << "better minimum found! check intial conditions and fitting parameters!" << std::endl; 
            new_min->PrintResults();
            exit(1);
          }
          if (new_min->MinValue() > min_chi2 + 1) {
            last_gt_val = val;
            last_gt_chi2 = new_min->MinValue();
            //interpolate
            double valtmp = last_lt_val + (last_gt_val - last_lt_val) * ((min_chi2 + 1 - last_lt_chi2)/(last_gt_chi2 - last_lt_chi2));
            if (std::abs(valtmp - last_gt_val) / std::abs(valtmp - last_lt_val) < 0.2 ) { 
              val = valtmp - 0.1*std::abs(valtmp - last_lt_val); 
            }
            else { 
              val = valtmp; 
            }
            //val = val - stepsize;
            if (val < LL.at(i)) { val = LL.at(i); zeroct += 1; }
            else if (val > UL.at(i)) { val = UL.at(i); zeroct += 1; }
            else { zeroct = 0; }
            //if (last_dir > 0) { stepsize = stepsize/2; }
            //last_dir = -1;
          }
          else {
            if (last_gt_chi2 < 0) {
              val = val + stepsize;
            }
            else {
              last_lt_val = val;
              last_lt_chi2 = new_min->MinValue();
              //interpolate
              double valtmp = last_lt_val + (last_gt_val - last_lt_val) * ((min_chi2 + 1 - last_lt_chi2)/(last_gt_chi2 - last_lt_chi2));
              if (std::abs(valtmp - last_lt_val) / std::abs(valtmp - last_gt_val) < 0.2 ) { 
                val = valtmp + 0.1*std::abs(valtmp - last_gt_val); 
              }
              else { 
                val = valtmp; 
              }
            } 
            if (val < LL.at(i)) { val = LL.at(i); zeroct += 1; }
            else if (val > UL.at(i)) { val = UL.at(i); zeroct += 1; }
            else { zeroct = 0; }
            //if (last_dir < 0) { stepsize = stepsize/2; }
            //last_dir = 1;
          }
          if (zeroct >= 5) { break; }
          ct -= 1;
        } 
        high = val - opt_val;
        //lower error
        if (opt_val < 1e-5) { low = -opt_val; }
        else {
          std::cout << "finding lower limit" << std::endl;
          last_gt_chi2 = -1;
          last_lt_chi2 = min_chi2;
          last_gt_val = -1;
          last_lt_val = opt_val;
          val = (opt_val + min_val)/2.;
          stepsize = (opt_val - min_val)/2.;
          last_dir = 0;
          ct = 100;
          zeroct = 0; 
          while (ct > 0) {
            new_min->SetFixedVariable(i, names[i], val);
            new_min->Minimize();
            //min->PrintResults();
            std::cout << std::endl << opt_val << "  " << val << "   " << min_val << "  " << max_val << "  " << new_min->MinValue() << "  " << min_chi2 << std::endl;
            if (std::abs(new_min->MinValue() - (min_chi2 + 1)) < 0.001) { break; }
            else if (new_min->MinValue() < min_chi2 - 1.0) { 
              std::cout << "better minimum found! check intial conditions and fitting parameters!" << std::endl; 
              new_min->PrintResults();
              exit(1);
            }
            if (new_min->MinValue() > min_chi2 + 1) {
              last_gt_val = val;
              last_gt_chi2 = new_min->MinValue();
              //interpolate
              double valtmp = last_lt_val + (last_gt_val - last_lt_val) * ((min_chi2 + 1 - last_lt_chi2)/(last_gt_chi2 - last_lt_chi2));
              if (std::abs(valtmp - last_gt_val) / std::abs(valtmp - last_lt_val) < 0.2 ) {
                val = valtmp + 0.1*std::abs(valtmp - last_lt_val);
              }
              else { 
                val = valtmp; 
              }
              //val = val + stepsize;
              //if (last_dir < 0) { stepsize = stepsize/2; }
              //last_dir = 1;
              if (val < LL.at(i)) { val = LL.at(i); zeroct += 1; }
              else if (val > UL.at(i)) { val = UL.at(i); zeroct += 1; }
              else { zeroct = 0; }
            }
            else {
              if (last_gt_chi2 < 0) {
                val = val - stepsize;
              }
              else {
                last_lt_val = val;
                last_lt_chi2 = new_min->MinValue();
                //interpolate
                double valtmp = last_lt_val + (last_gt_val - last_lt_val) * ((min_chi2 + 1 - last_lt_chi2)/(last_gt_chi2 - last_lt_chi2));
                if (std::abs(valtmp - last_lt_val) / std::abs(valtmp - last_gt_val) < 0.2 ) {
                  val = valtmp - 0.1*std::abs(valtmp - last_lt_val);
                }
                else { 
                  val = valtmp; 
                }
              }
              //val = val - stepsize;
              if (val < LL.at(i)) { val = LL.at(i); zeroct += 1; }
              else if (val > UL.at(i)) { val = UL.at(i); zeroct += 1; }
              else { zeroct = 0; }
              //if (last_dir > 0) { stepsize = stepsize/2; }
              //last_dir = -1;
            }
            if (zeroct >= 5) { break; }
            ct -= 1;
          } 
          low = val - opt_val; 
          std::cout << std::endl << "found my own errors of " << opt_val << " " << low << " " << high << std::endl;
         }
        delete new_min;
      } 
      min->SetPrintLevel(0);

      
      ampfile << i << "   " << names.at(i) << "   " << min->X()[i] << "   " << low << "    " << high << std::endl;
      if (i > 7 ) {
        ff->low_errors[i-8] = low;
        ff->up_errors[i-8] = high;
      }
    }
    else if (!nofit) {
      ampfile << i << "   " << names.at(i) << "   " << min->X()[i] << "   " << min->Errors()[i] << std::endl;
      if (i > 7 ) {
        ff->low_errors[i-8] = min->Errors()[i];
        ff->up_errors[i-8] = min->Errors()[i];
      }
    }

  }
  ampfile.close();
   
  TGraph *gr = ff->getfunc(min->X());
  std::ofstream datfile(outdatfile.c_str());
  datfile << "#fitted output from fitGretina" << std::endl;
  for (int i=0; i<gr->GetN(); ++i) {
    double x, y;
    gr->GetPoint(i, x,y);
    datfile << x << "  " << y << std::endl;
  }
  datfile.close();

  std::cout << Form("Writing histogram and graphs to %s", outrootfile.c_str()) << std::endl;
  TFile *outfile = new TFile(outrootfile.c_str(), "recreate");

  TH1D *expthist_bs = ff->get_bs_hist();
  expthist_bs->SetDirectory(outfile);
  expthist_bs->Write();
  TCanvas *c1 = new TCanvas();
  std::vector<int> colors={kGreen+2, kBlue, kOrange+1, kMagenta+1, kBlack, kRed+2, kGreen, kGray};
  expthist_bs->Draw("hist");
  expthist_bs->GetXaxis()->SetRangeUser(low, high);
  gr->SetLineWidth(2);
  gr->Draw("same L");
  gr->SetName("fitTotal");
  gr->Write();
  gr->SetEditable(false);    
  for (int i=0; i<ff->scales.size(); ++i) {
    //std::cout << ff->names[i] << "   " << ff->scales[i] << std::endl;
    if (!ff->types[i].compare("NUC")) {
      TGraph *gr_i = ff->getfunc(i, min->X());
      gr_i->SetName(ff->names[i].c_str());
      gr_i->SetLineColor(colors[i%8]);
      gr_i->Draw("same L");
      gr_i->SetEditable(false);
      gr_i->Write();
    }
  }
  TGraph *gr_bkg = ff->getbgfunc(min->X());
  gr_bkg->SetName("neutron_background");
  gr_bkg->SetLineColor(kGreen+2);
  gr_bkg->Draw("same L");
  gr_bkg->SetEditable(false);
  gr_bkg->Write();



  c1->Print((outpdf+"(").c_str());

  expthist_bs->Draw("E1 same");
  gr->Draw("same L");

  std::cout << "Evaluating chi^2 in ROI (" << ROI_ll << ", " << ROI_ul << ")" << std::endl;
  
  expthist_bs->GetXaxis()->SetRangeUser(ROI_ll, ROI_ul);
  //evaluate chi2 in ROI
  double chi2 = 0;
  int ndf = 0;
  double maxcont = 0;
  if (!method.compare("LOGLIKELIHOOD")) {
    double nll = 0;
    for (int i=expthist_bs->GetXaxis()->FindBin(ROI_ll); i<=expthist_bs->GetXaxis()->FindBin(ROI_ul); ++i) {
      double x, y;
      gr->GetPoint(i-1,x,y);
      if (expthist_bs->GetBinContent(i) > maxcont) {
        maxcont = expthist_bs->GetBinContent(i);
      }
      y+= ff->np_expthist->GetBinContent(i)*ff->np_scale;
      double yi = ff->expthist->GetBinContent(i);
      if (yi > 0) {
        y = std::max(y, 1e-5);
        nll += (y + yi * TMath::Log(yi/y) - yi);
        //std::cout << x << "  " << ff->expthist->GetXaxis()->GetBinCenter(i) << "  " << y << "  " << yi << "   " << (y + yi * TMath::Log(yi/y) - yi) << std::endl;
      }
      else {
        nll += y;
      }
    
      ndf += 1;
    }
    chi2 = 2.0*nll;
  }
  else if (!method.compare("CHI2")){
    for (int i=expthist_bs->GetXaxis()->FindBin(ROI_ll); i<=expthist_bs->GetXaxis()->FindBin(ROI_ul); ++i) {
      if (expthist_bs->GetBinContent(i) > maxcont) {
        maxcont = expthist_bs->GetBinContent(i);
      }
      std::pair<double, double> val;
      val = ff->eval(expthist_bs->GetBinCenter(i));
      if (expthist_bs->GetBinContent(i) > 0) {
        chi2 += std::pow((expthist_bs->GetBinContent(i) - val.first), 2)/(std::pow(expthist_bs->GetBinError(i), 2) + std::pow(val.second, 2));
      }
      else {
        chi2 += std::pow((expthist_bs->GetBinContent(i) - val.first), 2)/(1.0 + std::pow(val.second, 2));
      }
      ndf += 1;
    }
  }
  std::cout << "Chi^2 in ROI = " << chi2 << ", npoints = " << ndf << std::endl;
  
  expthist_bs->GetYaxis()->SetRangeUser(0, maxcont*1.2);
  TLatex *text = new TLatex();
  text->SetText(0.15, 0.85, Form("#chi^{2} = %4.1f, #chi^{2}/#nu = %4.2f", chi2, chi2/ndf));
  text->SetNDC();
  text->Draw();
    
  c1->Print(outpdf.c_str());

  TGraph *nondop_gr = ff->nondop_getfunc(min->X());  

  TH1D *expthist_nd = ff->get_nd_hist();
  expthist_nd->SetDirectory(outfile);
  expthist_nd->Write();
  TCanvas *c2 = new TCanvas();
  ff->nondop_expthist->Draw("hist");
  ff->nondop_expthist->GetXaxis()->SetRangeUser(low, high);
  nondop_gr->SetLineWidth(2);
  nondop_gr->Draw("same L");
  nondop_gr->SetEditable(false);
  nondop_gr->SetName("nondopFitTotal");
  nondop_gr->Write();
  for (int i=0; i<ff->scales.size(); ++i) {
    //std::cout << ff->names[i] << "   " << ff->scales[i] << std::endl;
    TGraph *gr_i = ff->nondop_getfunc(i, min->X());
    gr_i->SetName(("nondop_"+ff->names[i]).c_str());
    gr_i->SetLineColor(colors[i%8]);
    gr_i->Draw("same L");
    gr_i->SetEditable(false);
    gr_i->Write();
  }
  c2->Print((outpdf+")").c_str());
    
  //return min->MinValue();
  outfile->Close();
  delete outfile;
  return chi2;
}  

