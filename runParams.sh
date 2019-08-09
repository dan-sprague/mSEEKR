#!/bin/bash

for a in .99 .999 .9999
do
  for b in .99 .999 .9999
  do
    for c in noDog noRat noPig noHuman noMouse noBoar noMarm noCow
    do
      python train.py --query ./counts/${c}AE.skr --null ./counts/hg26Transcripts.skr -k 2,3,4,5 --qPrefix ${c}AE --nullPrefix hg26Trscpt --qT $a --nT $b --dir ./markovModels/${a:1}_${b:1}/
    done
  done
done
