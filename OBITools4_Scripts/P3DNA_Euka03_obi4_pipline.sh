#!/bin/bash
#SBATCH --job-name=P3_DNA_Euka03_obi4
#SBATCH --partition=cpu
#SBATCH --output=/work/FAC/FBM/DEE/lfumagal/default/julie/thesis/scripts/STD/P3_DNA_Euka03_obi4_output
#SBATCH --error=/work/FAC/FBM/DEE/lfumagal/default/julie/thesis/scripts/STD/P3_DNA_Euka03_obi4_error
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=julie.guenat@unil.ch
obi_dir="/users/jguenat/obi4/bin"
export PATH="$obi_dir:$PATH"
# 1. Sequences alignment
obipairing  --min-identity=0.8 --min-overlap=10 -F /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/data/P3DNA_F.fastq -R /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/data/P3DNA_R.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/1.obi4_illupairend_P3DNA_Euka03.fastq
#2. remove unaligned sequences
obigrep -p 'annotations.mode!="join"' /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/1.obi4_illupairend_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/2.obi4_rm_unalignedseq_P3DNA_Euka03.fastq
#3. Assign seq to samples
obimultiplex -t /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/data/P3DNA_Euka03_ngsfilter.txt -u /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/obi4_unidentified_P3DNA_Euka03.fastq /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/2.obi4_rm_unalignedseq_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/3.obi4_ngsfilter_assigned_P3DNA_Euka03.fastq
#4. dereplicate reads into uniq sequences
#grouping striclty identical reads together
obiuniq -m sample  /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/3.obi4_ngsfilter_assigned_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/4.obi4_dereplicated_seq_P3DNA_Euka03.fastq
#5. keep only key arguments
obiannotate -k count -k merged_sample /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/4.obi4_dereplicated_seq_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/5.obi4_keyargu_P3DNA_Euka03.fastq
#6. Tag seq for PCR errors
obiclean -s sample -r 0.05 -H /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/5.obi4_keyargu_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/6.obi4_cleanPCR_P3DNA_Euka03.fastq
#7. Keeping sequences having a count greater or equal to 10
obigrep -l 49 -L 264 -p 'sequence.Count()>=10' /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/6.obi4_cleanPCR_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/7.obi4_count10_size_P3DNA_Euka03.fastq
#8.Convert into Obi1 to run sumaclust
obiconvert --output-OBI-header /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/7.obi4_count10_size_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/8.obi1_count10_size_P3DNA_Euka03.fastq
#9.SUMACLUST
/users/jguenat/sumaclust_v1.0.36/sumaclust /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/8.obi1_count10_size_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/9.obi1_sumaclust_P3DNA_Euka03.fastq
#10. Convert into obi4 again
obiconvert  --output-json-header /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/9.obi1_sumaclust_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/10.obi4_sumaclust_P3DNA_Euka03.fastq
#11. Taxonomic assignation
obitag -t /work/FAC/FBM/DEE/lfumagal/default/julie/databases/TAXO -R /work/FAC/FBM/DEE/lfumagal/default/julie/databases/20240213/Euka03/DB.Euka03.3.40.270.fasta /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/10.obi4_sumaclust_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/11.Taxo_assigned_P3DNA_Euka03.fastq
#12. adding the taxonomical path
obiannotate -t /work/FAC/FBM/DEE/lfumagal/default/julie/databases/TAXO --taxonomic-path /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/11.Taxo_assigned_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/12.Taxo_pathadded_P3DNA_Euka03.fastq
#13. convert into obi1
obiconvert --output-OBI-header /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/12.Taxo_pathadded_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/13.Obi1_Taxo_pathadded_P3DNA_Euka03.fastq
#14.convert it into a tab file
source /users/jguenat/mambaforge/etc/profile.d/conda.sh
conda activate OBI
obitab -o /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/13.Obi1_Taxo_pathadded_P3DNA_Euka03.fastq > /work/FAC/FBM/DEE/lfumagal/default/julie/thesis/results/P3DNA_Euka03/14.tab_P3DNA_Euka03.txt
