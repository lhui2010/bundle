
Step By Step
1. Remove duplicated genes with samtools faidx
2. Remove strange characters like : or |
3. Add suffix to genes and chromosomes
4. Clean additional characters that follows space in fasta header
5. Translate cds to prot by replacing U to X

Exceptions:
1. Some GFF do not use CDS to represent genes
2. Sometimes gffread fails
