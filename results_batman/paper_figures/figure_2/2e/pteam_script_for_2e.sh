for i in {1..10..1}
do

for j in 1 2 3 4 5 7 9
do

for k in {0..4..1}
do

qsub -cwd -l m_mem_free=1G -b y Rscript --vanilla run_pTEAM_for_2e.R $i $j $k

done
done
done
