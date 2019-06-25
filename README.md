# mSEEKR: Markov based Sequence Evaluation through k-mer representation

## Example usage
python mSEEKR.py  --model ./markovModels/ --db ../markovAligner/mamXists.fa --nRAND 1000 --prefix mamxist -p .01 -w 200 -s 20 -n 7 --bkg ../standard_sequences/gencode.vM17.lncRNA_transcripts.fa

Parameter explanation:<br/>
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
