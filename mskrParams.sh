#!/bin/bash
for w in 50 100 200 500
do
  for k in 2 3 4 5
  do
    for c in AE BCD
    do
      python mSEEKR.py --model ./markovModels/mam${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/okRsx.fa --prefix rsx -w $w
    done
  done
done
