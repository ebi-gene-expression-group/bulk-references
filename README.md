# bulk-references
Collection of scripts to manage references for bulk-rnaseq. It uses nf-core/references to manage references and genome indices for tools like HISAT and STAR.

## Parabricks STAR references

Set `PARABRICKS_STAR_INDEX=true` to apply the Parabricks STAR override used by
the RNA-seq workflow. Its STAR genome-generation process then uses STAR
`2.7.2a`, the version required by the Parabricks `rna_fq2bam` container.
