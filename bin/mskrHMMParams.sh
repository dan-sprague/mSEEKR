#!/bin/bash

for a in 99 999 9999
do
  for b in 99 999 9999
  do
    for c in A B D E
    do
      python mSEEKR.py --db ./fastaFiles/mamX/human.fa --model ./markovModels/${a}_${b}_${c}_hgT/ -k 2,3,4,5,6,7,8,9 --prefix hgX
     done
  done
done
