#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input  = params.input  ?: 'star_inputs.tsv'
params.output = params.output ?: 'star_profiles.config'

def readRows(String path) {
    def lines = new File(path).readLines().findAll { it?.trim() }
    assert lines.size() > 1 : "Input TSV is empty or missing data: ${path}"

    def header = lines[0].split('\t', -1) as List
    lines.drop(1).collect { line ->
        def vals = line.split('\t', -1) as List
        def row = [:]
        header.eachWithIndex { h, i -> row[h] = i < vals.size() ? vals[i] : '' }
        row
    }
}

def asLong(row, key, long defaultValue = 0L) {
    def v = row[key]?.toString()?.trim()
    return v ? v.toLong() : defaultValue
}

def sanitizeProfileName(String s) {
    s.toLowerCase()
     .replaceAll(/[^a-z0-9]+/, '_')
     .replaceAll(/^_+|_+$/, '')
}

def suggestMultimap(long genomeSizeBp, long nContigs, long n50Bp) {
    int mm
    if (genomeSizeBp < 20_000_000L) {
        mm = 10
    } else if (genomeSizeBp < 200_000_000L) {
        mm = 20
    } else if (genomeSizeBp < 1_000_000_000L) {
        mm = 50
    } else {
        mm = 100
    }

    // Fragmented assemblies often benefit from a more permissive cap
    if (n50Bp < 100_000L || nContigs > 10_000L) {
        mm = Math.max(mm, 50)
    }
    return mm
}

def roundIntronMax(long proposed) {
    if (proposed <= 5_000L) return 5_000
    if (proposed <= 10_000L) return 10_000
    if (proposed <= 20_000L) return 20_000
    if (proposed <= 50_000L) return 50_000
    if (proposed <= 100_000L) return 100_000
    if (proposed <= 200_000L) return 200_000
    if (proposed <= 500_000L) return 500_000
    return 1_000_000
}

def suggestAlignIntronMax(long maxIntronBp, long p99IntronBp) {
    if (maxIntronBp <= 0L && p99IntronBp <= 0L) {
        return 10_000
    }
    long proposed = Math.max((long)(p99IntronBp * 1.5d), (long)(maxIntronBp * 1.1d))
    return roundIntronMax(proposed)
}

def suggestMotifs(long genomeSizeBp, long maxIntronBp) {
    // Conservative rule of thumb:
    // compact/simple genomes with short introns -> RemoveNoncanonical
    // larger/complex/non-model genomes -> None
    if (genomeSizeBp < 50_000_000L && maxIntronBp > 0L && maxIntronBp < 10_000L) {
        return 'RemoveNoncanonical'
    }
    return 'None'
}

workflow {
    def rows = readRows(params.input)

    def config = new StringBuilder()

    config << "params {\n"
    config << "  star_base_args = [\n"
    config << "    '--quantMode TranscriptomeSAM',\n"
    config << "    '--outSAMtype BAM Unsorted',\n"
    config << "    '--outSAMattributes NH HI AS NM MD',\n"
    config << "    '--readFilesCommand zcat',\n"
    config << "    '--twopassMode Basic',\n"
    config << "    '--runRNGseed 0',\n"
    config << "    '--alignSJDBoverhangMin 1',\n"
    config << "    '--outSAMstrandField intronMotif',\n"
    config << "    '--quantTranscriptomeSAMoutput BanSingleEnd'\n"
    config << "  ].join(' ')\n"
    config << "}\n\n"

    config << "profiles {\n"

    rows.each { row ->
        String dataset = row['dataset'] ?: 'dataset'
        String profile = sanitizeProfileName(dataset)

        long genomeSizeBp = asLong(row, 'genome_size_bp')
        long nContigs     = asLong(row, 'n_contigs')
        long n50Bp        = asLong(row, 'n50_bp')
        long maxIntronBp  = asLong(row, 'max_intron_bp')
        long p99IntronBp  = asLong(row, 'p99_intron_bp')

        int mm            = suggestMultimap(genomeSizeBp, nContigs, n50Bp)
        int intronMax     = suggestAlignIntronMax(maxIntronBp, p99IntronBp)
        String motifs     = suggestMotifs(genomeSizeBp, maxIntronBp)

        config << "  ${profile} {\n"
        config << "    process {\n"
        config << "      withName: 'STAR_ALIGN' {\n"
        config << "        ext.args = [\n"
        config << "          params.star_base_args,\n"
        config << "          '--outFilterMultimapNmax ${mm}',\n"
        config << "          '--alignIntronMax ${intronMax}',\n"
        config << "          '--outFilterIntronMotifs ${motifs}'\n"
        config << "        ].join(' ')\n"
        config << "      }\n"
        config << "    }\n"
        config << "  }\n\n"
    }

    config << "}\n"

    def out = new File(params.output)
    out.text = config.toString()

    println "Wrote ${out.absolutePath}"
}
