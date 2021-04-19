#!/bin/sh

# DESCRIPTION
#   compare if the gfortran and Python version
#   yield the same result

basedir=`pwd`
GMAP_fortran_url="https://github.com/IAEA-NDS/GMAP-Fortran.git"
GMAP_fortran_commit_id_01="78c669b9fd6fdf4b57e73bb89030134a1b323bb1"
GMAP_fortran_commit_id_02="43b42145b53d2379f686f488f7c27eb906945c14"

run_fortran_gmap() {

    GMAP_fortran_commit_id=$1
    GMAP_idx=$2

    if [ ! -d "GMAP-Fortran/$GMAP_idx" ]; then
        mkdir -p "GMAP-Fortran/$GMAP_idx"
        git clone $GMAP_fortran_url GMAP-Fortran/$GMAP_idx
        retcode=$?
        if [ $retcode -ne "0" ]; then
            print "Downloading Fortran version of GMAP from $GMAP_fortran_url failed"
            exit $retcode
        fi
        if [ ! -d "GMAP-Fortran "]; then
            print "The directory expected to contain the Fortran source code does not exist"
            exit $?
        fi
        cd GMAP-Fortran/$GMAP_idx \
        && git checkout $GMAP_fortran_commit_id \
        && rm -rf .git \
        && cd source \
        && gfortran -o GMAP GMAP.FOR \
        && mv GMAP .. \
        && chmod +x GMAP
        cd "$basedir"
    fi

    if [ ! -d result/fortran ]; then
        mkdir -p result/$GMAP_idx
        cp input/data.gma result/$GMAP_idx/
        cd result/$GMAP_idx/
        ../../GMAP-Fortran/$GMAP_idx/GMAP
        # Fortran print -0.0 whereas Python 0.0
        # thus remove the odd minus sign from Fortran output
        # two times each sed comment to deal with the case
        # -0.00-0.00 and to make sure both minuses are removed
        sed -i -e 's/-\(00*\.0*\)\([^0-9]\|$\)/ \1\2/g' gma.res
        sed -i -e 's/-\(00*\.0*\)\([^0-9]\|$\)/ \1\2/g' gma.res
        sed -i -e 's/-\(00*\.0*\)\([^0-9]\|$\)/ \1\2/g' plot.dta
        sed -i -e 's/-\(00*\.0*\)\([^0-9]\|$\)/ \1\2/g' plot.dta
        cd "$basedir"
    fi
}



# generate Python result
run_fortran_gmap $GMAP_fortran_commit_id_01 01
run_fortran_gmap $GMAP_fortran_commit_id_02 02

# compare the results
cd "$basedir"
diff result/01/gma.res result/02/gma.res \
  && diff result/01/plot.dta result/02/plot.dta

retcode=$?
exit $retcode
