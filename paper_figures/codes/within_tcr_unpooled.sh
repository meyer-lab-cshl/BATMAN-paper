for i in {1..70..1}
do
for j in {1..66..1}
do

qsub -cwd -l m_mem_free=1G -b y Rscript --vanilla within_tcr_unpooled.R $i $j

done
done
