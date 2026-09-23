# bulk-references
Collection of scripts to manage references for bulk-rnaseq. It uses nf-core/references to manage references and genome indices for tools like HISAT and STAR.

## Parabricks STAR references

Set `PARABRICKS_STAR_INDEX=true` to apply the Parabricks STAR override used by
the RNA-seq workflow. Its STAR genome-generation process then uses STAR
`2.7.2a`, the version required by the Parabricks `rna_fq2bam` container, and
writes the versioned STAR index alongside the other reference outputs. With
the variable unset, the historical reference-generation behavior is retained.

Run it from a Codon environment with the reference variables configured:

```bash
PARABRICKS_STAR_INDEX=true sbatch create_reference_with_indices.sh datasheet.yaml
```

The generated Parabricks index is expected under the versioned STAR reference
directory ending in `2.7.2a`. Keep the existing `2.7.11b` index for native STAR
workflows; the two indices must not be mixed.
