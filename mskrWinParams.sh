#!/bin/bash

for w in 50 100 200 500
do
  for k in 2 3 4 5
  do
    for c in AE BCD
    do
      python mSEEKR.py --model ./markovModels/mam${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/oRsx.fa --prefix orsx -w $w
      python mSEEKR.py --model ./markovModels/mam${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/kRsx.fa --prefix krsx -w $w
    done
  done
done
