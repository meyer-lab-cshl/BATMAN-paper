for i in {1..100..1}
do

export PYTENSOR_FLAGS="compiledir=$HOME/.pytensor/compiledir_7_1_$i"
qsub -cwd -l m_mem_free=1G -b y python loo_tcr_AL_variable_steps_npeptide.py $i 7 1


done
