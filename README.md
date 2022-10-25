# ASGAL-DROS

```
conda create -c bioconda -n psiesg suppa biopython gffutils pysam samtools gffread shark rich
conda activate psiesg
git clone --recursive https://github.com/ldenti/asgaldros.git
cd asgaldros/galig
git submodule update --init --recursive
git checkout snake+psi
make prerequisites
make
cd ..
```

```
cd example
tar xvfz data.tar.gz
python3 ../psiesg.py -l 7 --threads 4 reference.fa annotation.gtf sample.fq
```

### TODO
* [X] ~extend to A3 and A5~
* [X] ~code refactoring~
* [X] ~extend to IR and other events from SUPPA2~
* [X] ~extend to MX~
* [ ] allow user to compute PSI from .bam
* [ ] what if constitutive exons are shorter than read length?
* [ ] improve PSI (e.g., use exon coverage)
* [ ] preserve gene strand
