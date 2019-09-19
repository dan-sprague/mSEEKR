mSEEKR

##################################
Installation
##################################
  1. Dependencies
    a. Python3.6
      check with command "python -V"

      if < 3.6, easiest solution to download anaconda python3.6/3.7
    b. cython
      Try "which cython", if none

      pip install cython

      OR

      conda install -c anaconda cython

  2. Type following commands in desired directory

    git clone https://github.com/spragud2/mSEEKR.git
    cd mSEEKR/
    python setup.py build_ext --inplace

  3. Ignore warnings (unless further steps don't work)

##################################
Counting k-mers
##################################

  1. Curate unique fasta files for queries and null model before hand
  2. Use the following command

  python kmers.py --fasta ./fastaFiles/mouseA.fa -k 2,3,4 --name mouseA -n 3
  python kmers.py --fasta ./fastaFiles/gencode.vM17.lncRNA_transcripts.fa -k 2,3,4 --name mm10Trscpts -n 3


  Parameters:

    --fasta : path to fasta file
    -k : comma separated list of values of k
    --name : name to give output file
    -n : number of processing cores to use (scales with value of k -- so should be less than or equal to number of k values specified)


  Output:

  Outputs binary .skr (seekr) files containing count matrices to /counts/ directory

##################################
Training markov models
##################################

  0. Count k-mers for queries and null before proceeding
  1. Run the following command

  python train.py --query ./counts/mouseA.skr --null ./counts/mm10Trscpts.skr -k 2,3,4 --qPrefix mouseA --nullPrefix mm10Trscpts --qT .9999 --nT .9999 --dir ./markovModels/

  Parameters:

  --query : Path to query count file
  --null : Path to null count file
  -k : comma separated list of values of k to train for, must have been calculated prior
  --qPrefix : prefix file name for query
  --nullPrefix : prefix file name for null

  --qT : Query to query transition parameter, query to null is 1 - qT
  --nT : Null to null transition parameter, null to query is 1 - nT

  --dir : output directory

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

###################################################
Find HMM state path through sequences of interest
###################################################
  1. Run the following command

  python mSEEKR.py --db ./fastaFiles/mamX/mouse.fa -n 8 --prefix test --model ./markovModels/


  Parameters

  --db : sequences of interest to run the HMM on
  --model : path to directory containing models (ie /markovModels/)
  -n : Number of processing cores to use. Scales with size of fasta file (# of entries, not sequence length)
  --prefix : experiment prefix name
