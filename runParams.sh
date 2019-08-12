#!/bin/bash

for c in noDog noRat noPig noHuman noMouse noBoar noMarm noCow mam
do
  for d in AE BCD
  do
    python train.py --query ./counts/${c}${d}.skr --null ./counts/hg26Transcripts.skr -k 2,3,4,5,6 --qPrefix ${c}${d} --nullPrefix hg26Trscpt
  done
done
