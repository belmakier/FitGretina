#include "FitGretina.hh"

int main(int argc, const char **argv) {
  if (argc < 2) {
    std::cout << "useage: fitGretina [input file] [output file]" << std::endl;
    return 0;
  }
  fitManager fm;
  std::string inpfile(argv[1]);
  fm.Load(inpfile);
  double minchi2 = fm.DoFit();
  std::ofstream outfile(argv[2], std::ios::app);
  outfile << "   " << minchi2 << std::endl;
  outfile.close();
    
}
