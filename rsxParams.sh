#!/bin/bash

for c in mamAE mamBCD
do
  python train.py --query ./counts/${c}.skr --null ./counts/hg26Transcripts.skr -k 2,3,4,5 --qPrefix ${c} --nullPrefix hg26Trscpt
done
