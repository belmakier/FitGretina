#include "FitGretina.hh"

int main(int argc, const char **argv) {
  fitManager fm;
  std::string inpfile(argv[1]);
  fm.Load(inpfile);

  std::vector<double> ppar_lo;
  std::vector<double> ppar_hi;

  double lo = std::atof(argv[2]);
  double hi = std::atof(argv[3]);
  int npoints = std::atoi(argv[4]);
  int thispoint = -1;
  if (argc>5) {
    thispoint = std::atoi(argv[5]); 
  }

  for (int i=0; i<npoints; ++i) {
    ppar_lo.push_back(lo + (double)i*(hi-lo)/(double)npoints);
    ppar_hi.push_back(lo + (double)(i+1)*(hi-lo)/(double)npoints);
  }

  std::vector<std::vector<double> > amplitudes;  //indexed by state, point
  std::vector<std::vector<double> > low_errors;
  std::vector<std::vector<double> > up_errors;
  for (int i=0; i<fm.types.size(); ++i) {
    if (fm.types[i].compare("NUC")) { continue; }
    std::vector<double> v;
    amplitudes.push_back(v);
    low_errors.push_back(v);
    up_errors.push_back(v);
  }

  for (int i=0; i<npoints; ++i) {
    if (thispoint>0 && i<thispoint) {
      for (int j=0; j<fm.ff->scales.size(); ++j) {
        if (fm.ff->types[j].compare("NUC")) { continue; }
        amplitudes[j].push_back(0);
        low_errors[j].push_back(0);
        up_errors[j].push_back(0);
      }        
    }
    else {
      fm.gate_low = ppar_lo[i];
      fm.gate_high = ppar_hi[i];
      fm.Setup();
      fm.DoFit();
      for (int j=0; j<fm.ff->scales.size(); ++j) {
        if (fm.ff->types[j].compare("NUC")) { continue; }
        amplitudes[j].push_back(fm.ff->scales[j]);
        low_errors[j].push_back(fm.ff->low_errors[j]);
        up_errors[j].push_back(fm.ff->up_errors[j]);
      }        
    }
  }

  for (int i=0; i<amplitudes.size(); ++i) {
    std::ofstream distrfile(Form("ppar_distr_%s.dat", fm.ff->names[i].c_str()));
    for (int j=0; j<npoints; ++j) {
      if (thispoint>0 && j<thispoint) { continue;}
      distrfile << j << "   " << ppar_lo[j] << "   " << ppar_hi[j] << "   " << amplitudes[i][j] << "   " << low_errors[i][j] << "   " << up_errors[i][j] << std::endl;
    }
    distrfile.close();
  }
}

    
    
