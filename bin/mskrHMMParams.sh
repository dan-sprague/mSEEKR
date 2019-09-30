#!/bin/bash

for a in 99 999 9999
do
  for b in 99 999 9999
  do
    for c in A B D E
    do
      #python mSEEKR.py --model ./markovModels/${a}_${b}/${c}_hg26Trscpt --db ./fastaFiles/hgkcnq1ot1.fa --prefix ${a}_${b}_kcn
      python mSEEKR.py --db ./fastaFiles/mamX/human.fa --model ./markovModels/${a}_${b}_${c}_hgT/ -k 2,3,4,5,6 --prefix hgX
    #    python mSEEKR.py --model ./markovModels/${a}_${b}/noHuman${c}_hg26Trscpt --db ./fastaFiles/mamX/human.fa --prefix ${a}_${b}
    #    python mSEEKR.py --model ./markovModels/${a}_${b}/noRat${c}_hg26Trscpt --db ./fastaFiles/mamX/rat.fa --prefix ${a}_${b}
    #    python mSEEKR.py --model ./markovModels/${a}_${b}/noPig${c}_hg26Trscpt --db ./fastaFiles/mamX/pig.fa --prefix ${a}_${b}
    #    python mSEEKR.py --model ./markovModels/${a}_${b}/noMouse${c}_hg26Trscpt --db ./fastaFiles/mamX/mouse.fa --prefix ${a}_${b}
    #    python mSEEKR.py --model ./markovModels/${a}_${b}/noDog${c}_hg26Trscpt --db ./fastaFiles/mamX/dog.fa --prefix ${a}_${b}
    #    python mSEEKR.py --model ./markovModels/${a}_${b}/noCow${c}_hg26Trscpt --db ./fastaFiles/mamX/cow.fa --prefix ${a}_${b}
    #    python mSEEKR.py --model ./markovModels/${a}_${b}/noMarm${c}_hg26Trscpt --db ./fastaFiles/mamX/marm.fa --prefix ${a}_${b}
     done
  done
done
