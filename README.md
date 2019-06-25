# mSEEKR: Markov based Sequence Evaluation through k-mer representation

## Workflow

1. Train transition matrices using train.py
2. Do this for any number of models desired (different sequences, different values of -k, etc...)
3. If you desire statistical evaluation, run mSEEKR.py 
4. If you desire sequence plots, run makePlots.py


## Example usage

### train.py 
python train.py --query ./fastaFiles/mamxist_BCD_lnc_4_200w_20sl.fa --null ./fastaFiles/gencode.vM17.lncRNA_transcripts.fa --qPrefix mamBCD4 --nullPrefix lnc -k 3,4

#### Parameter explanation:<br/>
--query,type=str,help='Path to fasta file or containing sequences to build markov model (e.g. functional regions of a ncRNA)<br/>
--null, type=str,help='Path to fasta file containing sequences that compose null model (e.g. transcriptome, genome, etc.)<br/>
--qPrefix,type=str,help='String, Output file prefix',default=query<br/>
--nullPrefix,type=str,help='String, Output file prefix',default=null<br/>
--dir,type=str,help='Output directory for matrix files',default=./markovModels/<br/>
-k,type=str,help='Comma delimited string of possible k-mer values',default=2,3,4<br/>
-a,type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default=ATCG<br/>

### mSEEKR.py 
python mSEEKR.py  --model ./markovModels/ --db ../markovAligner/mamXists.fa --nRAND 1000 --prefix mamxist -p .01 -w 200 -s 20 -n 7 --bkg ../standard_sequences/gencode.vM17.lncRNA_transcripts.fa

#### Parameter explanation:<br/>
--model: Path to directory containing .mkv files OR path to single .mkv file<br/>
--db: Path to fasta file containing sequences to compare to<br/>
--nRAND: Number of random sequences to generate for statistical model<br/>
--prefix: Add this prefix to output files<br/>
--retrain: Save unique sequence hits from running this program against --db to a fasta file.<br/> 
--bkg: Path to fasta file to calculate background nucleotide frequencies. Default is uniform<br/>
-p: 0 < p < 1, Desired p-val of log-likelihood ratio score "S" to consider significant. If P(S > 0) is very small and less than this argument, set S = 0 and p = P(S>0); default=.01<br/>
-k: k-mer length<br/>
-w: window length<br/>
-s: slide size<br/>
-a: Alphabet. Default = ATCG<br/>
-n: number of CPU cores to use (useful for large databases)<br/>

### makePlots.py 

python makePlots.py --db ../markovAligner/mamXists.fa -e mamXists --win 100,200,300,400,500


#### Parameter explanation:<br/>
--mkv,type=str,help='Directory containing markov (.mkv) models',default=./markovModels/<br/>
--db,type=str,help='Path to fasta file with sequences to calculate similarity score'<br/>
--win,type=str,help='Comma delimited string (no quotes) of window sizes,default=100,200<br/>
--slide,type=str,help='Comma delimited string (no quotes) of how much to slide for each window size specified. If a single number is specified, that slide will be applied to all windows, otherwise, the length of this list must equal the window argument length',default=20<br/>
--dir,type=str,help='Path to save directory,default=./pdfs/<br/>
-a,type=str,help='Alphabet, string (no quotes)',default=ATCG<br/>
-e,type=str,help='Experiment name (no quotes)'<br/>
