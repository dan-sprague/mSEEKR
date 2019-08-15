#!/bin/bash

for w in 50 100 200 500
do
  for k in 2 3 4 5 6
  do
    for c in A B D E
    do
      python mSEEKR.py --model ./markovModels/mam${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/oRsx.fa --prefix orsx -w $w
      python mSEEKR.py --model ./markovModels/mam${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/kRsx.fa --prefix krsx -w $w
      python mSEEKR.py --model ./markovModels/noHuman${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/mamX/human.fa --prefix human -w $w
      python mSEEKR.py --model ./markovModels/noMouse${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/mamX/mouse.fa --prefix mouse -w $w
      python mSEEKR.py --model ./markovModels/noMarm${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/mamX/marm.fa --prefix marmoset -w $w
      python mSEEKR.py --model ./markovModels/noDog${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/mamX/dog.fa --prefix dog -w $w
      python mSEEKR.py --model ./markovModels/noCow${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/mamX/cow.fa --prefix cow -w $w
      python mSEEKR.py --model ./markovModels/noRat${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/mamX/rat.fa --prefix rat -w $w
      python mSEEKR.py --model ./markovModels/noPig${c}_hg26Trscpt_${k}.mkv --db ./fastaFiles/mamX/pig.fa --prefix pig -w $w
    done
  done
done
