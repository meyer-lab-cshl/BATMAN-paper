for i in {0..151..1}
do


export PYTENSOR_FLAGS="compiledir=$HOME/.pytensor/compiledir_${i}"
qsub -cwd -l m_mem_free=1G -b y python run_BATMAN_for_S5a.py $i

done

