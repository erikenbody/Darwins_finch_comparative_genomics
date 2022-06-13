
#for every sample, run one job per chromosome for whatshap
#script found in the tool directory 

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_whatshap
cd $WORK_D

while read SAMPLE
do
  echo "$SAMPLE"
  mkdir -p $SAMPLE
  cd $SAMPLE

  sbatch --array=1-31 ~/afc/10_LDhelmet/tools/submit_parv_whatshap_V2.sh $SAMPLE

  cd $WORK_D

done < parv_names.txt
