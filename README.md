# ASGAL-DROS

```
create -c bioconda -n asgaldros suppa biopython gffutils pysam samtools gffread shark snakemake-minimal

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

### TODO
* [ ] what if constitutive exons are shorter than 100bp?
* [ ] extend to A3, A5, IR and other events from SUPPA2