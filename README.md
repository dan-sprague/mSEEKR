# mSEEKR


## Installation

### Dependencies

#### Python3.6
check with command 

```
python -V
```

if < 3.6, easiest solution is to download anaconda python3.6/3.7
#### cython

Try 
```
which cython
``` 
if none

```
pip install cython
```

OR

```
conda install -c anaconda cython
```

#### Type following commands in desired directory
```
	git clone https://github.com/spragud2/mSEEKR.git
	cd mSEEKR/
	python setup.py build_ext --inplace
```
#### Ignore warnings (unless further steps don't work)
<hr/>

## Counting k-mers 


  1. Curate unique fasta files for queries and null model before hand
  2. Use the following command
```
  python kmers.py --fasta ./fastaFiles/mouseA.fa -k 2,3,4 --name mouseA -n 3
  python kmers.py --fasta ./fastaFiles/gencode.vM17.lncRNA_transcripts.fa -k 2,3,4 --name mm10Trscpts -n 3
```

#### Parameters:

<ol>
1. --fasta : path to fasta file
2. -k : comma separated list of values of k
3. --name : name to give output file
4. -n : number of processing cores to use (scales with value of k -- so should be less than or equal to number of k values specified)
5. --dir : output directory 


  Output:

  Outputs binary .skr (seekr) files containing count matrices to --dir

<hr/>
## Training markov models

  0. Count k-mers for queries and null before proceeding
  1. Run the following command
```
  python train.py --query ./counts/mouseA.skr --null ./counts/mm10Trscpts.skr -k 2,3,4 --qPrefix mouseA --nullPrefix mm10Trscpts --qT .9999 --nT .9999 --dir ./markovModels/
```

Parameters:

1. --query : Path to query count file
2. --null : Path to null count file
3. -k : comma separated list of values of k to train for, must have been calculated prior
4. --qPrefix : prefix file name for query
5. --nullPrefix : prefix file name for null
6. --qT : Query to query transition parameter, query to null is 1 - qT
7. --nT : Null to null transition parameter, null to query is 1 - nT
8. --dir : output directory

  Output:

  Directory containing models for each value of k specified, the directory structure would look like:

    | markovModels
    |
    |--- mouseA_mouseNull
    |------- 2
    |------------- model.hmm
    |------------- model.mkv
    |------- 3
    |------------- model.hmm
    |------------- model.mkv
    |------- 4
    |------------- model.hmm
    |------------- model.mkv
    |--- mouseB_mouseNull
    .
    .
    .

```
Find HMM state path through sequences of interest
```
  1. Run the following command

  python mSEEKR.py --db ./fastaFiles/mamX/mouse.fa -n 8 --prefix test --model ./markovModels/


  Parameters

  1. --db : sequences of interest to run the HMM on
  2. --model : path to directory containing models (ie /markovModels/)
  3. -k : value of k to use in the analysis (must have been specified in training)
  4. -n : Number of processing cores to use. Scales with size of fasta file (# of entries, not sequence length)
  5. --prefix : experiment prefix name

