#!/bin/bash
for a in .5 .6 .7 .8 .9 .99 .999 .9999
do
  for b in .5 .6 .7 .8 .9 .99 .999 .9999
  do
    for c in A B D E
      do
      python train.py --query ./counts/${c}.skr --null ./counts/hgT.skr -k 2,3,4,5,6 --qPrefix ${a:1}_${b:1}_${c} --nPrefix hgT --qT ${a} --nT ${b} --dir ./markovModels/
    done
  done
done


#noPigA noPigB noPigD noPigE noCowA noCowB noCowD noCowE noRatA noRatB noRatD noRatE noDogA noDogB noDogD noDogE noMarmA noMarmB noMarmD noMarmE noMouseA noMouseB noMouseD noMouseE noHumanA noHumanB noHumanD noHumanE
