# mSEEKR

This is a program for identifying regions of high similarity based on *k*-mer content to some set of query sequences, relative to a null background set of sequences.

## Dependencies

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
## Installation

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
  python kmers.py --fasta ./fastaFiles/mA.fa -k 2,3,4 --name mouseA -n 3 --dir ./counts/
  python kmers.py --fasta ./fastaFiles/gencode.vM17.lncRNA_transcripts.fa -k 2,3,4 --name mm10Trscpts -n 3 --dir ./counts/
```

#### Parameters:

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
  mv *skr counts/
  python train.py --query ./counts/mouseA.skr --null ./counts/mm10Trscpts.skr -k 2,3,4 --qPrefix mouseA --nPrefix mm10Trscpts --qT .9999 --nT .9999 --dir ./markovModels/
```

Parameters:

1. --query : Path to query count file
2. --null : Path to null count file
3. -k : comma separated list of values of k to train for, must have been calculated prior
4. --qPrefix : prefix file name for query
5. --nPrefix : prefix file name for null
6. --qT : Query to query transition parameter, query to null is 1 - qT
7. --nT : Null to null transition parameter, null to query is 1 - nT
8. --dir : output directory

  Output:

  Directory containing models for each value of k specified, the directory structure would look like:

    | markovModels
    |
    |--- mouseA_mouseNull
    |------- 2
    |------------- hmm.mkv
    |------- 3
    |------------- hmm.mkv
    |------- 4
    |------------- hmm.mkv
    |--- mouseB_mouseNull
    .
    .
    .

## (OPTIONAL) Optimize transition parameters using the B-W Algorithm
This script will take the output .mkv output file from train.py and find an MLE for the transition parameters. 

  1. Run the following command
```
  python bw.py -k 4 --db ./fastaFiles/xist.fa --prior markovModels/mouseA_mm10Trscpts/4/hmm.mkv --its 3
```

The above command takes the HMM trained on repeat A at k =4, and uses Xist to find an MLE for the transition parameters. More than one sequence can be provided as training, if a multi-entry FASTA file is provided. 

Parameters:

1. -k : value of k for k-mers
2. --db : path to FASTA file to optimize parameters on
3. --prior : hmm.mkv file output by train.py
4. --its : number of iterations of the baum-welch algorithm to run


Output:
1. Replaces the file passed as --prior with an updated version that contains the MLE transition matrix
2. A file containing the trajectory of the transition parameters over each loop of the BW algorithm, check this for convergence




## Find HMM state path through sequences of interest

  1. Run the following command
```
  python mSEEKR.py --db ./fastaFiles/mm10kcn.fa -n 8 --prefix test --model ./markovModels/mouseA_mm10Trscpts -k 3 --fasta
```

  Parameters

  1. --db : sequences of interest to run the HMM on
  2. --model : path to directory containing models (ie /markovModels/query_null/)
  3. -k : value of k to use in the analysis (must have been specified in training)
  4. -n : Number of processing cores to use. Scales with size of fasta file (# of entries, not sequence length)
  5. --prefix : experiment prefix name
  
