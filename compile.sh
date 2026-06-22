#!/bin/bash
set -e

DIR=$(realpath $(dirname $0))
OUTDIR=$(pwd)

stuff=$1

compstr=""
compstr_FCOV=""
if [[ "$stuff" == *"g"* ]] ; then
	if [[ "$stuff" == *"p"* ]] ; then
		compstr="-O0 -ggdb -fsanitize=address -fsanitize=undefined -DDEBUG -ffree-line-length-400 -finit-local-zero -finit-real=zero -finit-integer=0"
	else
		compstr="-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow -ggdb -DDEBUG -ffree-line-length-400"
	fi
	compstr_FCOV=$compstr
else
	compstr="-O3 -funroll-loops -march=native -ffree-line-length-400"
	compstr_FCOV="-O2 -funroll-loops -fopenmp-simd -fopt-info-vec-missed -march=native -ffree-line-length-400"
fi

if [[ "$stuff" == *"w"* ]] ; then
	compstr="${compstr} -fmax-errors=3"
fi

ompstr=""
if [[ "$stuff" == *"o"* ]] ; then
	ompstr="-fopenmp"
fi

profstr=""
#if [[ "$stuff" == *"p"* ]] ; then
#	if [[ "$stuff" == *"r"* ]] ; then	
#		profstr="-pg -O0 -g"
#	elif [[ "$stuff" == *"c"* ]] ; then
#		profstr="-O0 -g"
#	else
#		profstr="-p -O0 -g"
#	fi
#fi

faststr=""
fastt="no"
if [[ "$stuff" == *"f"* ]] ; then
	faststr="-ffast-math"
	fastt="yes"
fi

#Does not work on GCC
#Might be that I did not set it up correctly or GCC does not support both OpenMP and OpenACC
if [[ "$stuff" == *"a"* ]] ; then
	accstr="-B /storage/praha1/home/zdrahalj/my_packages/gcc-12-nvptx/usr/lib/gcc/x86_64-linux-gnu/12/ -fopenacc -fopt-info-optimized-omp -foffload=nvptx-none -foffload-options=nvptx-none=-lfopenacc"
fi


echo "DIR"
echo $DIR
echo "OUTDIR"
echo $OUTDIR
echo "compiled with:"
echo $compstr
echo "openmp:"
echo $ompstr
echo "profiler:"
echo "openacc:"
echo $accstr
echo $profstr
echo "fast math? ${fastt}"

set -v
gfortran -w ${compstr} -c -cpp $DIR/constants.f95 -o $OUTDIR/constants.o
gfortran -w ${compstr} -c -cpp $DIR/strings.f95 -o $OUTDIR/strings.o
gfortran -w ${compstr} -c -cpp $DIR/util.f95 -o $OUTDIR/util.o
gfortran -w ${compstr} -cpp $DIR/suben.f95 -o $OUTDIR/suben
gfortran -w ${compstr} -c -cpp $DIR/wrram.f95 -o $OUTDIR/wrram.o

if [[ "$stuff" != *"j"* ]] ; then
	f95 -w -O3 $DIR/gar9.f -o $OUTDIR/gar9
	f95 -w -O3 $DIR/new1.f -o $OUTDIR/new1
	f95 -w ${compstr} $DIR/new2.f -o $OUTDIR/new2
	f95 -w ${compstr} $DIR/new4.f -o $OUTDIR/new4
	gfortran -w ${compstr} $DIR/geten_ZPE.f95 -ffree-line-length-200 -o $OUTDIR/geten_ZPE
	gfortran -w ${compstr} $OUTDIR/util.o $OUTDIR/strings.o $OUTDIR/constants.o $DIR/dusch_rev1.f $DIR/lapack_dgetrf.f -o $OUTDIR/dusch_rev1
	gfortran -w ${compstr} -cpp $DIR/grtoq.f95 -o $OUTDIR/grtoq
	gfortran -w ${compstr} $OUTDIR/constants.o $DIR/geten.f95 -ffree-line-length-200 -o $OUTDIR/geten
	gfortran -w ${compstr} -cpp $DIR/pmz_rev1.f -o $OUTDIR/pmz_rev1
	gfortran -w ${compstr} $DIR/gar_orca.f95 -o $OUTDIR/gar_orca
	gfortran -w ${compstr} -cpp $OUTDIR/strings.o $DIR/pmz_orca.f95 -o $OUTDIR/pmz_orca
	gfortran -w ${compstr} -cpp $OUTDIR/util.o $OUTDIR/strings.o $OUTDIR/constants.o $DIR/tms_sumRules.f95 -o $OUTDIR/tms_sumRules

	gfortran -w ${compstr} -cpp -ffree-line-length-500 $OUTDIR/constants.o $OUTDIR/wrram.o $OUTDIR/strings.o $OUTDIR/util.o $DIR/FCOV_spectrum.f95 -o $OUTDIR/FCOV_spectrum
	gfortran -w ${compstr} -cpp $OUTDIR/constants.o $OUTDIR/strings.o $OUTDIR/util.o $DIR/tdd2_num.f95 -o $OUTDIR/tdd2_num
	gfortran -w ${compstr} ${ompstr} -cpp $OUTDIR/constants.o $OUTDIR/strings.o $OUTDIR/util.o $OUTDIR/wrram.o $DIR/rroa_td_num.f95 -o $OUTDIR/rroa_td_num
	gfortran -w ${compstr} -cpp $OUTDIR/constants.o $OUTDIR/strings.o $DIR/gettms.f95 -o $OUTDIR/gettms
	gfortran -w ${compstr} -cpp $OUTDIR/constants.o $OUTDIR/strings.o $OUTDIR/util.o $DIR/elpol.f95 -o $OUTDIR/elpol
	gfortran -w ${compstr} -cpp $OUTDIR/constants.o $OUTDIR/util.o $OUTDIR/strings.o $DIR/addpol.f95 -o $OUTDIR/addpol
fi

gcc -c ${compstr} $DIR/cpu.c -o $OUTDIR/cpu.o
gfortran -cpp -c ${compstr} ${ompstr} ${accstr} ${profstr} ${faststr} -ffree-line-length-500 $DIR/big_numbers.f95 -o ${OUTDIR}/big_numbers.o
rm -rf lapack_temp
mkdir lapack_temp
cd lapack_temp
gfortran -cpp -c -O2 $DIR/lapack/util/*.f $DIR/lapack/lapack_routine/*.f
cd ..
gfortran -cpp -fmax-errors=10 ${compstr_FCOV} ${ompstr} ${accstr} ${profstr} ${faststr} -fopenmp-simd -fopt-info-vec -ffree-line-length-500 $OUTDIR/constants.o $OUTDIR/strings.o $OUTDIR/cpu.o $OUTDIR/util.o $DIR/FCOV_new.f95 $OUTDIR/lapack_temp/*.o ${OUTDIR}/big_numbers.o -o $OUTDIR/FCOV

rm -r lapack_temp
rm $OUTDIR/*.o
rm *.mod
set +v
exit



