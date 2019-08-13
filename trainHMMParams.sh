#!/bin/bash
for a in .99 .999 .9999
do
  for b in .99 .999 .9999
  do
    for c in humanAE humanBCD
    do
      python train.py --query ./counts/${c}.skr --null ./counts/hg26Transcripts.skr -k 2,3,4,5,6 --qPrefix ${c} --nullPrefix hg26Trscpt --qT $a --nT $b --dir ./markovModels/${a:1}_${b:1}/
    done
  done
done
