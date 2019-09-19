#!/bin/bash
for a in .99 .999 .9999
do
  for b in .99 .999 .9999
  do
    for c in hgD
      do
      python train.py --query ./counts/${c}.skr --null ./counts/hgT.skr -k 2,3,4,5,6 --qPrefix ${c} --nullPrefix hgT --qT $a --nT $b --dir ./markovModels/${a:1}_${b:1}
    done
  done
done


#noPigA noPigB noPigD noPigE noCowA noCowB noCowD noCowE noRatA noRatB noRatD noRatE noDogA noDogB noDogD noDogE noMarmA noMarmB noMarmD noMarmE noMouseA noMouseB noMouseD noMouseE noHumanA noHumanB noHumanD noHumanE
