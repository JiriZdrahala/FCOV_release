# Programs for calculation of resonance Raman and ROA

## Vibronic calculation
The vibronic calculation is done via the single-state approach. That means the electronic states and their vibrations are treated one at a time and negligible nonadiabatic coupling is present in the system. Time dependent (TD) approach is the main focus now, time independent (TI) is only kept for reference and validation. \
2nd derivatives of transition dipole/quadrupole moments can be included in the calculation, although they need to be calculated numerically from 1st derivatives so far.
In the TD approach, for a molecule with 93 vibrational modes, 2^16 time grid points, on 8 CPU threads the runtime of the calculation of the cross-correlation function is:\
**1st derivatives**: 1 minute\
**2nd derivatives**: 5 minutes, 38 seconds (not including the calculation of the derivatives themselves)\
\
In the .docx files there lies the derivation of the time-dependent approach. It is methodology heavy and unpolished, but after reading it (and with the help of LLM) you should be able to replicate or expand the TD approach. If you see any errors/oversights, please message me or post an issue.

## Installation
Go to your favourite programs folder and call the script "compile.sh o". If you want debugging, compile with "compile.sh og"

## Index of notable programs
FCOV_new.f95 - SINGLE-STATE vibronic calculation of resonance Raman and ROA with the TI and TD approach\
FCOV_spectrum.f95 - creation of RR and RROA spectra from .POLARS files\
addpol.f95 - summing polarizabilities in .POLARS files\
dusch_rev1.f - Duschinsky transformation, originally written by prof. Petr Bouř, revisioned by Jiří Zdráhala\
gar9, new1, new2, new4.f - information extraction from Gaussian output files, written by prof. Petr Bouř\
pmz_rev1.f - genereation of displaced geometries for numerical differentiation, originally written by prof. Petr Bouř, revisioned by Jiří Zdráhala\
rroa_td_num.f95 - calculation of resonance Raman and ROA through numerical differentiation of polarizabilities.\
tdd2_num.f95 - calculation of second (and third diagonal) derivatives of transition moments, done numerically\

## Where do I find...?
*...cross-correlation function calculation?* - FCOV_new.f95, subroutine Make_Corrf_SplitPropagator(...)\
\
*...Franck Condon integral calculation?* - FCOV_new.f95, function MakeFCArr_recalculate(...)\
*...excited vibrational state prescreening?* - FCOV_new.f95, function DEM_1v_RROA(...)

## Actually useful sources and not just lazy mentions
Lami, A. and Santoro, F. (2011). Time-Dependent Approaches to Calculation of Steady-State Vibronic Spectra: From Fully Quantum to Classical Approaches. In Computational Strategies for Spectroscopy, V. Barone (Ed.). https://doi.org/10.1002/9781118008720.ch10

Alberto Baiardi, Julien Bloino, Vincenzo Barone; General Time Dependent Approach to Vibronic Spectroscopy Including Franck–Condon, Herzberg–Teller, and Duschinsky Effects. J. Chem. Theory Comput. 10 September 2013; 9 (9): 4097–4115. https://doi.org/10.1021/ct400450k

Alberto Baiardi, Julien Bloino, Vincenzo Barone; A general time-dependent route to Resonance-Raman spectroscopy including Franck-Condon, Herzberg-Teller and Duschinsky effects. J. Chem. Phys. 21 September 2014; 141 (11): 114108. https://doi.org/10.1063/1.4895534

Alberto Baiardi, Julien Bloino, Vincenzo Barone; Time-Dependent Formulation of Resonance Raman Optical Activity Spectroscopy. J. Chem. Theory Comput. 11 December 2018; 14 (12): 6370–6390. https://doi.org/10.1021/acs.jctc.8b00488

Bernardo de Souza, Giliandro Farias, Frank Neese, Róbert Izsák; Efficient simulation of overtones and combination bands in resonant Raman spectra. J. Chem. Phys. 7 June 2019; 150 (21): 214102. https://doi.org/10.1063/1.5099247
