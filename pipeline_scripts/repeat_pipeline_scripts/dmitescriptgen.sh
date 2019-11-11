#!/bin/bash

for i in "$@"
do
case $i in
    -d=*|--dir=*)
    GENDIR="${i#*=}"
    shift # past argument=value
    ;;
    -n=*|--name=*)
    GENNAME="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--prefix=*)
    PREFIX="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--cpus=*)
    CPUS="${i#*=}"
    shift # past argument=value
    ;;
    -s=*|--settings=*)
    SETTINGS="${i#*=}"
    shift # past argument=value
    ;;
    -dm=*|--dmdir=*)
    DMDIR="${i#*=}"
    shift # past argument=value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
          # unknown option
    ;;
esac
done
echo "GENOME DIRECTORY  = ${GENDIR}"
echo "GENOME NAME       = ${GENNAME}"
echo "PROJECT PREFIX    = ${PREFIX}"
echo "CPUS              = ${CPUS}"
echo "NAME OF SETTINGS  = ${SETTINGS}"
echo "DETECT MITE DIR   = ${DMDIR}"
echo "FULL PATH TO GENOME = ${GENDIR}/${GENNAME}"
if [[ -n $1 ]]; then
    echo "Extra unnecessary argument supplied. This was ignored."
fi
# Generate output file
echo "clear;clc;" > ${SETTINGS}
echo "data_file = '${GENDIR}/${GENNAME}';" >> ${SETTINGS}
echo "genome_name = '${PREFIX}';" >> ${SETTINGS}
echo "cd ${DMDIR}" >> ${SETTINGS}
echo >> ${SETTINGS}
echo "tic;" >> ${SETTINGS}
echo "        do_MITE_detection(data_file,'-genome',genome_name,'-cpu',${CPUS})" >> ${SETTINGS}
echo "runtime = toc;" >> ${SETTINGS}
echo >> ${SETTINGS}
echo "fid = fopen('detectMITE.Runtime.txt','wt');" >> ${SETTINGS}
echo "fprintf(fid,'-------%s-------\n',genome_name);" >> ${SETTINGS}
echo "fprintf(fid,'Runtime: %f s\n',runtime);" >> ${SETTINGS}
echo "fclose(fid);" >> ${SETTINGS}
echo "exit;" >> ${SETTINGS}
