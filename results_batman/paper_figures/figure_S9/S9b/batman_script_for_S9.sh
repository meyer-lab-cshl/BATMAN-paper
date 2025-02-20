for i in {1..50..1}
do
for j in {0..11..1}
do

export PYTENSOR_FLAGS="compiledir=$HOME/.pytensor/compiledir_${i}_${j}"
qsub -cwd -l m_mem_free=1G -b y python run_BATMAN_AL_random_after_2_steps_for_S9.py $i $j

done
done
