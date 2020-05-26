#!/bin/sh

# DESCRIPTION
#   compare if the gfortran and ifort compiled 
#   executables yield the same result

gfortran_exe=gfortran
ifort_exe=ifort

basedir=`pwd`
fortran_src="$basedir/../../source/GMAP.FOR"

rm -r $basedir/result/gfortran
rm -r $basedir/result/ifort
mkdir -p $basedir/result/gfortran
mkdir -p $basedir/result/ifort
cp $basedir/input/data.gma $basedir/result/gfortran/
cp $basedir/input/data.gma $basedir/result/ifort/

echo
echo ----- Compiling and running the gfortran version of GMAP -----
echo
cd $basedir/result/gfortran
$gfortran_exe -o GMAP_gfortran $fortran_src \
&& ./GMAP_gfortran > gmap_gfortran.out \

echo
echo ----- Compiling and running the ifort version of GMAP -----
echo
cd $basedir/result/ifort
$ifort_exe -o GMAP_ifort $fortran_src \
&& ./GMAP_ifort > gmap_ifort.out \

cd $basedir/result
diff ifort/gma.res gfortran/gma.res > diff_result.out
ret=$?

if [ "$ret" -ne "0" ]; then
    echo
    echo Differences between gfortran and ifort version detected.
    echo Check result/diff_result.out for the specific differences.
fi

exit $ret
