scp -r Desktop/aziz/diseaseMechanism/* mezlinim@scc1.bu.edu:/project/lipidext/diseaseMechanism/

1-Preparation:
-Install the followig R packages: Matrix, preprocessCore
source("http://bioconductor.org/biocLite.R"); biocLite("preprocessCore");
-Create a data directory under /project/lipidext/diseaseMechanism/ where you will put your data and annotations. you will need these 5 files into it: 

1- A patient annotation file ("patients.txt")
2- A variants annotation file ("variants.txt")
3- An Exome input file ("genotype.txt")
4- A gene expression input file ("gene_expression.txt")
5- An optional gene annotation file ("genes.txt") 


-Starting from an Annovar annotated VCF file, we provide a script (vcf_to_inputs.py) that generate the genotype.txt and the variants.txt files with the correct format.
The VCF should be first annotated with Annovar with the following command:
#table_annovar.pl -protocol refGene,ljb26_all,snp138,esp6500siv2_all -operation g,f,f,f -vcfinput -buildver hg19
Then the script can be used as follows:
vcf_to_inputs.py $annovar_multianno_vcf_file $output_directory

2-Patients annotation:
Prepare a file listing all patients and their phenotypes. 2 columns: patients names(ids), phenotypes (0=controls , 1=cases)
The phenotypes can also be continuous between 0 and 1 (indicating the probability of being a case). Individuals that are not very close to 0 or 1 will have much weaker effect on the analysis. 
We recommend that the data should be balanced, i.e with as many cases as controls, or a symetric distribution around 0.5 in the continuous case.
Any patients whose names are not in this file will be ignored in exome and expression data. 


3-Variants Annotation:

We provide a script (vcf_to_inputs.py) to automatically generate this file from an Annovar-annotated VCF file.
Only the coding variants should go into this file: nonsense (splicing, frameshift, stop gain/loss) and missense, and they all should have associated harmfulness predictions.
The annotation file should have these tab-separated columns:
variant_id, Gene, Type, Harmfulness, OtherColumns
Where variant id is of the format "chr_start_Ref_Alt" (chromosome number, start position, reference allele, alternative allele), type is for example in [splicing, stop gain, stop loss, frameshift, nonsynonymous], harmfulness are predictions between 0(neutral) and 1(damaging) obtained from Polyphen2 for example. Nonsense variants should be given a high harmfulness (0.9 by default).
Other columns can be added at the end: for example here we added MAF, which can be missing for some variants (blank). 
Variants ids that are not in this file will be ignored in exome data.

Here is an example of annotation file:

"10_52587953_C_A"	"A1CF"	nonsynonymous	0.79	0.11
"10_52595854_G_A"	"A1CF"	frameshift	0.95	
"10_52596055_C_T"	"A1CF"	nonsynonymous	0.75	0.06
"10_52601632_A_T"	"A1CF"	nonsynonymous	0.35	0.01
"12_9220358_-_T"	"A2M"	stop-gain	0.95	
"12_9221429_G_A"	"A2M"	nonsynonymous	0.95	0.15
"12_9242994_C_T"	"A2M"	nonsynonymous	0.70	0.22


4-Exome Input:
We provide a script (vcf_to_inputs.py) to automatically generate this file from an Annovar-annotated VCF file.
This file contains all heterozygous and homozygous alternative calls, in all patients. The file should contain three columns: Variant_id (same as in variants annotation), patient_id (same as in patients annotation), zygosity (1=heterozygous, 2=homozygous Alt). 
Only include the variants in the annotation file (nonsense and missense) as everything else will be ignored. Reference calls (genotype =0) should not be included.
 
N.B: Be careful and verify that there is no discrepancies between the variant ids in this file and in the variants.txt file. Sometimes Annovar changes the REF/ALT values, especially for indels. This is important since we need reliable keys to map variants in the exome inputs files with their annotations in the annotation file (chr_start_Ref_Alt is a good key if Ref and Alt unchanged).


5-Gene Expression inputs:
Prepare a file containing the matrix of gene expression levels. The patients ids are the column names and gene ids are the rown names. 
We use the raw expression data. If preprocessed data is available (corrected for batch effects or quantile normalized) you should skip that operation from our code so that preprocessing is not done twice.
The file should be directly readable by R into the correct data.frame and the gene names should generally match with the ones in the genes/variants annotation and the patients ids should match with the ones in the patients annotation.

6-Gene annotation:
Optional file containing a list of genes (one gene name per line). The analysis will be performed on these given genes, and will ignore all other genes present in the exome data, expression data or gene network. A second colum with priors over the genes can be added. If this file is not given, the list of genes present in the variants annotation file will be used by default.

7-Run the method:
The file analyse_data.r can be used directly to analyse your method.
If the file naming convention is respected you can run your analysis with the following command:
analyse_data.r $diseaseMechanism_folder $data_folder $results_folder $exprfile $genefile"

where $exprfile and $genefile indicate whether you want to use expression data and whether a genes file is provided(1=TRUE (default filename), 0=FALSE, another_string=filename).

8-Examples:
We provide simulation data in the simul/simul_10_0_100/rep1/ directory. You will find the input files in the correct format and a file with the true simulated causal genes.
To test the method on that simulations set use the command: 
Rscript Rproject/analyse_data.r "./" "simul/simul_10_0_100/rep1/" "simul/simul_10_0_100/" 1 1

Example of a script to submit a job with analyse_data.r:
echo "/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/Rproject/analyse_data.r /hpf/largeprojects/agoldenb/aziz/diseaseMechanism/ /hpf/largeprojects/agoldenb/aziz/copd/ /hpf/largeprojects/agoldenb/aziz/copd/results/ 0 0" | qsub -l mem=63G,vmem=63G,nodes=1:ppn=1,walltime=23:50:00 -N "copd"
