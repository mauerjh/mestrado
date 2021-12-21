

while read p;
do /home/escience/anaconda3/bin/atropos --mirna -a AACTGTAGGCACCATCAAT -m 16 -M 35 -o trim_"$p" -se /genetica_2/Jessica/seq_teste/arquivos/"$p" -T 12 --report-file "$p".txt
done < /genetica_2/Jessica/seq_teste/arquivos/sample_file.txt

while read p;
do fastqc ./"$p"
done < ./list_files.txt

while read p;
do docker run -v /genetica_1/Jessica/mestrado/seqrenata_fastq/trim:/exceRptInput -v /genetica_1/Jessica/mestrado/seqrenata_fastq/align:/exceRptOutput -v /genetica_2/excerpt_files/hg38_dir:/exceRpt_DB/hg38 -t rkitchen/excerpt INPUT_FILE_PATH=/exceRptInput/"$p" ENDOGENOUS_LIB_PRIORITY=miRNA,tRNA,piRNA,circRNA,gencode N_THREADS=12 JAVA_RAM=32G ADAPTER_SEQ=none REMOVE_LARGE_INTERMEDIATE_FILES=true
done < ./list_files.txt