=================================================
                    FitGretina
 A library for fitting GRETINA spectra, in combination with UCGretina
 Tim Gray - tgray30@utk.edu
==================================================

Principles of operation
-----------------------
This library/program is designed to be used in the context of fast-beam experiments with GRETINA. The idea is that we populate some residual nucleus in a reaction, populating several excited states which decay by gamma emission. We want to understand the states that are populated - their properties (i.e., lifetimes, decay branches), as well as their relative populations.

To do this, we simulate (using UCGretina) a response function for each possible state in the residual nucleus. I direct you to the UCGretina github for how to do this. What we have at the end of this are

1. An experimental Doppler-corrected spectrum (or possibly multiple spectra, more on that below)
2. A collection of response functions (or histograms) representing what GRETINA will see for each of the states that could be populated.

What we want to do is appropriately scale the response function histograms so that their sum matches the experimental spectrum. At this point, the scaling factors represent the relative population amplitudes. This program is designed to assist with conducting this fit.

Background
----------

Unfortunately, the experimental spectrum is not only composed of response functions from the residual nucleus of interest. In reality, there are other things creating gamma rays that need to be accounted for. There are three primary sources which are accounted for in different ways by FitGretina

1. non-correlated background (i.e., room background)

These are gamma rays that happen to be detected by GRETINA but are not induced by the beam. Think K40 or gammas from natural radon decay series. These are accounted for by providing an experimental non-prompt spectrum: one where the gamma-ray and particle (residual in the S800 or light particle detected from another auxiliary detector) are not coincident --- they are separated by a relatively long time. FitGretina will add this non-prompt background, with a fixed scaling factor to the response functions.

2. neutron background

At least in the context of knockout reactions (for which this code was developed), a relatively large number of neutrons are produced by the beam. These then inelastically scatter and cause prompt gamma rays in reactions on various parts of the apparatus: Al and Fe in the target chamber/beam pipe etc. etc.

The gamma rays produced are prompt (i.e., correlated with the particle) and so will not be accounted for by the previous background. They can, however be simulated reasonably well by UCGretina, by providing a loosely forward-focused distribution of neutrons at various energies. 

The gamma rays themselves are generally not strongly Doppler-shifted: they are produced by nuclei that are relatively stationary in the lab frame. This means that in the experimental (i.e., Doppler-corrected) spectrum, they do not appear as distinct sharp transitions, but broad features that are not easily constrained. For this reason, FitGretina can fit both Doppler-corrected and non-Doppler-corrected spectra simultaneously: the non-Doppler-corrected spectrum shows sharp lines corresponding to the neutron background, and serves to constrain the magnitudes of this background.

Thus, for each response function/UCGretina simulation, we should provide both a Doppler-corrected and a non-Doppler-corrected spectrum. In addition, we need to provide both a Doppler-corrected and a non-Doppler-corrected experimental spectrum.

3. other background

In addition to the previous two backgrounds, an exponential function is often used as an empirical background. FitGretina allows for this, as well as the possibility for the sum of two exponential backgrounds with different decay constants. In my experience, the second exponential background is not needed if the neutron backgrounds are correctly accounted for.

Using FitGretina
----------------

FitGretina uses an input file to specify the fitting options, the experimental input, and the desired output. These input files are simple in format, with each line consisting of a keyword and then (usually) one or more values. The various keywords and their meaning are below

INPUT [filename]
The ROOT file which contains the experimental histograms.

TYPE [histogram type]
Type should be 1D or 2D corresponding to whether the experimental histograms are 1D or 2D. If 2D is used, one of two situations must be satisfied:
1. the PROJX keyword is used to specify which y-bin should be used. This is a useful mode if multiple residuals are being considered: i.e. the 2D histogram may be Doppler-corrected gamma energy (x) vs residual ID (y).
2. the GATE keyword is used to specify a y-range to gate on. This is useful when constructing, for example, a parallel-momentum distribution. A 2D histogram could be gamma energy (x) vs parallel momentum (y), and FitGretina could fit different slices of the parallel momentum distribution. 

DOP_P [histogram name]
Name of the prompt, Doppler-corrected histogram

DOP_NP [histogram_name]
Name of the non-prompt, Doppler-corrected histogram

NODOP [histogram_name]
Name of the non-Doppler-corrected histogram

PROJX [bin number]
Bin number to take an x-projection of.

GATE [low] [high]
Lower and upper bounds of the gate. Note that these correspond to y-axis units, not bins. The experimental spectra will be formed by projecting the 2D input histograms onto the x axis, after gating on the y range specified here.

NP_SCALE [scale]
Scaling factor for the non-prompt subtraction. Usually this will be the ratio of the prompt and non-prompt timing window widths.

ERROR [error type]
Specify how the fitted errors will be extracted. Possible values are MINOS or HESSE. HESSE is default, fast, but probably inaccurate. Good for getting a good fit dialed in, but don't trust the uncertainties. MINOS is more reliable and much slower. In reality, if MINOS is specified, I attempt to use the GetMinosError function from the ROOT minimizer, and if it fails I do my own estimation.

ROI [low] [high]
Region of interest for zoomed plot and chi^2. This is not used in the fit per se, but useful for inspecting the output.

FFILE [type] [filename] [histogram] [non-Doppler histogram] [label] [value]
This is the specification of a response function. 
type: possible values "NUC" or "BKG" corresponding to the residual of interest (NUC) or simulated backgrounds (BKG)
filename: filename in which the histogram(s) are found
histogram: Doppler-corrected response 
non-Doppler histogram: non-Doppler-corrected response
label: short, descriptive label to use for the minimizer
value: starting guess for parameter value

DOPBKG [amp1] [tau1] [amp2] [tau2]
starting parameters for the phenomenological exponential background. If amp2=tau2=0, then the second exponential will be turned off in the fit.

NDBKG [amp1] [tau1] [amp2] [tau2] 
Same as above but for the non-Doppler-corrected spectrum

RANGE [low] [high]
Fit range

OUTDAT [filename]
Output file (ASCII format) for final fitted curve

OUTROOT [filename]
Output ROOT file to create which will contain the experimental spectra as well as the fitted curve

AMPFILE [filename]
Output file for fitted amplitudes of each of the response functions

OUTPDF [filename]
Output file for summary PDF showing the Doppler-corrected and non-Doppler-corrected fits, as well as a zoomed portion for the ROI

REBIN [value]
Rebin factor for the experimental histograms

METHOD [type]
Valid methods are "CHI2" or "LOGLIKELIHOOD", which corresponds to the evaluation of the goodness of fit. For low statistics cases, make sure you use LOGLIKELIHOOD.

NOFIT
This has no keyword. If present, it will fix all parameters to their starting values and immediately return the fit (i.e., no adjustment or actual fitting takes place). This can be useful for picking good starting guesses, or adjusting things by hand to see the effect of raising or lowering particular states.

You can see some examples of input files in the examples/ folder

Fitting a single spectrum
-------------------------

The most basic task is fitting a single spectrum. A program which does this is contained with the library, it should be located in the bin/ directory. Its source can be viewed in src/SingleFit.cc

Its operation is very simple: it opens an input file specified by a command-line argument, and loads it into the "fitManager" object. The "fitManager" then conducts the fit with DoFit(), and writes the output final chi^2 to a file.

Fitting a momentum distribution
-------------------------------
Fitting a momentum distribution is a little more complex. An example of a program which does this can be found in src/pParFit.cc

Once again, we specify the input file with a command-line argument. This input file should include the GATE keyword, though we will specify the range ourselves in the pParFit program.  The (2D) histograms necessary for this should have gamma-ray energy on the x-axis, and then parallel momentum on the y-axis.

The second and third command line arguments correspond to the lower and upper bounds of the momentum distribution. The fourth corresponds to the number of momentum points or bins to use. The momentum range is split into equal bins, and the gate bounds in the fitManager object set manually. For each bin, a fit is conducted and the amplitudes and errors are stored. Finally, the amplitudes and errors are written out to text files for later plotting.

Hopefully these two programs give a sense of how a combination of the input file, command line arguments, and a custom program can be used flexibly to suit your particular fitting needs.

