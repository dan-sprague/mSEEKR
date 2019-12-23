#!/bin/bash

for a in 5 6 7 8 9 99 999 9999
do
  for b in 5 6 7 8 9 99 999 9999
  do
    for c in A B D E
    do
      python mSEEKR.py --db ./fastaFiles/hgkcnq1ot1.fa --model ./markovModels/${a}_${b}_${c}_hgT/ -k 2,3,4,5,6 --prefix hgKcn
     done
  done
done
