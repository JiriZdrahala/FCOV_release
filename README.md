# Programs for calculation of resonance Raman and ROA

## Vibronic calculation
2nd derivatives of transition dipole/quadrupole moments can be included in the calculation.
In the TD approach, for a molecule with 93 vibrational modes, 2^16 time grid points, on 8 CPU threads the runtime of the calculation of the cross-correlation function is:\
**1st derivatives**: 1 minute\
**2nd derivatives**: 5 minutes, 38 seconds (not including the calculation of the derivatives themselves)

## Installation
Go to your favourite programs folder and call "compile.sh o". If you want debugging, compile with "compile.sh og"

## Index of notable programs
FCOV_new.f95 - vibronic calculation of resonance Raman and ROA with the TI and TD approach\
FCOV_spectrum.f95 - creation of RR and RROA spectra from .POLARS files\
addpol.f95 - summing polarizabilities in .POLARS files\
dusch_rev1.f - Duschinsky transformation, originally written by prof. Petr Bouř, revisioned by Jiří Zdráhala\
gar9, new1, new2, new4.f - information extraction from Gaussian output files, written by prof. Petr Bouř\
pmz_rev1.f - genereation of displaced geometries for numerical differentiation, originally written by prof. Petr Bouř, revisioned by Jiří Zdráhala\
rroa_td_num.f95 - calculation of resonance Raman and ROA through numerical differentiation of polarizabilities.\
tdd2_num.f95 - calculation of second (and third diagonal) derivatives of transitiom moments, done numerically\

## Where do I find...?
*...cross-correlation function calculation?* - FCOV_new.f95, subroutine Make_Corrf_SplitPropagator(...)


## Useful sources
