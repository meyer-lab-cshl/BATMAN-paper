for i in {1..20..1}
do
for j in 1 2 3 4 5 7 9 10
do
for k in {0..4..1}
do

export PYTENSOR_FLAGS="compiledir=$HOME/.pytensor/compiledir_${i}_${j}_${k}"
qsub -cwd -l m_mem_free=1G -b y python run_unpooled_BATMAN_for_2e.py $i $j $k

done
done
done
