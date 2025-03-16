for i in {1..150..1}
do
for j in {0..11..1}
do

echo $i

rm -rf $HOME/.pytensor/compiledir_${i}_${j}

done
done

