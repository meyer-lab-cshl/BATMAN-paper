for i in {1..100..1}
do

export PYTENSOR_FLAGS="compiledir=$HOME/.pytensor/compiledir_$i"
qsub -cwd -l m_mem_free=1G -b y python loo_tcr_active_random_mixed.py $i

done
