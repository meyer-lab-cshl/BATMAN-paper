for i in {1..150..1}
do
for j in {0..11..1}
do

export PYTENSOR_FLAGS="compiledir=$HOME/.pytensor/compiledir_${i}_${j}"
qsub -cwd -l m_mem_free=1G -b y python run_BATMAN_AL_for_4a.py $i $j

done
done
