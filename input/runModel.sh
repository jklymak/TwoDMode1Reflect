#!/bin/sh -l
#PBS -m be
#PBS -M jklymak@gmail.com
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -l walltime=0:02:00
#PBS -q debug
#PBS -A ONRDC35552400
#PBS -j oe
#PBS -N ${JOBNAME}


# 14400 took 10 minutes, so this should get us to 288 h.  Probably not enough
# time, but there will be restarts available.  
#

cd $PBS_O_WORKDIR
# top=$1  Passed as qsub  -v top=h60h27f20 runModel.sh

PARENT=MarginalMixing
top=${PBS_JOBNAME}
results=${WORKDIR}/${PARENT}/
outdir=$results$top

# These should already be copied
#cp data00 $outdir/_Model/input/data
#cp eedata $outdir/_Model/input
#cp data.* $outdir/_Model/input
#cp ../build/mitgcmuv $outdir/_Model/build

#printf "Copying to archive"
#rm -rf ../archive/$top
#cp -r $outdir/_Model/ ../archive/$top
#rm -rf ../archive/$top/indata/*


pwd

rm STD*.*

if true; then
    ls -al mitgcmuv
    printf "Starting: $outdir\n"
    aprun -n 8 ./mitgcmuv > mit.out
fi


if false; then
# make mds files of last hour
cd $PBS_O_WORKDIR
echo ${top}
cd ../python
pwd
MDSJOB=`qsub -v PRE="${top}",ITER=720000,ITEREND=723600,DITER=1800 runmdstonetcdf.sh`

echo "Submited mds"
echo ${MDSJOB}

# now GetSlices:  This will wait for MDSJOB to finish...

cd $PBS_O_WORKDIR
cd ../python
pwd
qsub -v PRE="${top}",ITER=720000 -W depend=afterok:${MDSJOB} runGetSlice.sh

echo "Submited runGetSlice.sh"

# now do energy budget:  This will wait for MDSJOB to finish...

cd $PBS_O_WORKDIR
cd ../python
pwd
qsub -v PRE="${top}",U0="${U0}",ITER=720000,DITER=1800 -W depend=afterok:${MDSJOB} runGetEnergy.sh

echo "Submited getEnergy"

# Now archive
TODO=${PBS_JOBNAME}
cd ${WORKDIR}/${PARENT}
rm -f archive_job
cat > archive_job << END
#!/bin/bash
#PBS -m be
#PBS -M jklymak@gmail.com
#PBS -l walltime=24:00:00
#PBS -q transfer
#PBS -A ONRDC35552400
#PBS -l select=1:ncpus=1
#PBS -j oe
#PBS -S /bin/bash

echo "Transfer job ${PARENT}/${TODO} Started"
cd ${WORKDIR}/${PARENT}
tar cf ${TODO}_files.tar  ${WORKDIR}/${PARENT}/${TODO}
gzip ${TODO}_files.tar
rsh ${ARCHIVE_HOST} mkdir ${ARCHIVE_HOME}/${PARENT}
rcp ${TODO}_files.tar.gz  ${ARCHIVE_HOST}:${ARCHIVE_HOME}/${PARENT}
rsh ${ARCHIVE_HOST} ls -l  ${ARCHIVE_HOME}/${PARENT}
echo "Transfer job ${PARENT}/${TODO} ended"
END
#
# Submit the archive job script.
#qsub archive_job
# End of batch job


fi
