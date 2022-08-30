# ASGAL-DROS

```
conda create -c bioconda -n asgaldros suppa biopython gffutils pysam samtools gffread shark snakemake-minimal cutadapt pandas

wget http://ftp.flybase.net/releases/FB2022_04/dmel_r6.47/fasta/dmel-all-chromosome-r6.47.fasta.gz
wget http://ftp.flybase.net/releases/FB2022_04/dmel_r6.47/gtf/dmel-all-r6.47.gtf.gz

gunzip dmel-all-chromosome-r6.47.fasta.gz
samtools faidx dmel-all-chromosome-r6.47.fasta 2L 2R 3L 3R 4 X Y > dmel-all-chromosome-r6.47.chroms.fasta

gunzip dmel-all-r6.47.gtf.gz
grep -P "FBgn0000504|FBgn0003023|FBgn0264270|FBgn0003741|FBgn0003028" dmel-all-r6.47.gtf > sex-genes.gtf

# setup config.yaml
snakemake -s preprocessing.smk -j 4
snakemake -s asgal.smk -j 4
```

Compute PSI value from a BAM file
```
# run SUPPA2
suppa.py generateEvents -i [gtf] -o [prefix] -f ioe -e SE SS
# concatenate SUPPA2 .ioe
cat [prefix]_*.ioe > [all.ioe]
# compute PSI from ioe and bam
python3 compute_psi_bam.py [all.ioe] [bam]
```

### TODO
* [X] ~extend to A3 and A5~
* [ ] code refactoring
* [ ] extend to IR and other events from SUPPA2
* [ ] what if constitutive exons are shorter than read length?
* [ ] improve PSI
