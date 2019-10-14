#!/bin/bash
for a in 99 999 9999
do
  for b in 99 999 9999
  do
    for c in mamAE mamBCD
    do
      python mSEEKR.py --model ./markovModels/${a}_${b}/${c}_hg26Trscpt --db ./fastaFiles/okRsx.fa --prefix ${a}_${b}
    done
  done
done
