# Family tables `0.5.0`: what the tables do not show

`DESCRIBE` gives each table's columns and types. These notes hold what a query gets wrong without an error. A name
written with its table, as in `variants.genes.ensembl_gene_id`, is a column, or a field of a struct or of a list's
entries, and `<member>` stands for a member's column.

`manifest.json` lists the tables that the folder holds, the members (`sample_id`, `sequenced`, `affected`, `role`) and
the columns that the run left out because it did not name their source (`omitted_columns`, `sv_omitted_columns`). A
query on such a column fails.

## Rows

- `variants`: one row per alternate allele, after a site of several alleles is split. Rows count alleles, not sites.
- `transcripts`: one row per allele and transcript, of Ensembl and of RefSeq, so a join with `variants` multiplies rows.
  The Annotator flags MANE Select on the RefSeq transcript only: `transcripts.is_mane_select` with
  `transcripts.source = 'Ensembl'` matches nothing.
- `sv`: one row per record of two callers, told apart by `sv.caller` (`sv`, `cnv`). Both can report one event, and
  nothing links the two rows. A record that is also a row of `variants` names it in `sv.small_variant_id`, and `IS NULL`
  there keeps a query over both tables from counting it twice.
- `repeats`: one row per repeat of the catalogue, whatever the members' alleles. A row is not a finding.

## A member's column

`variants`, `sv` and `repeats` end in one struct column per sequenced member, named as the sample. Only the manifest
says who is the proband.

- `variants.<member>.gt` is the genotype after the split: `1/2` is `1/0` in one row and `0/1` in the other. Read
  `<member>.carries` and `<member>.zygosity` instead of matching it.
- `<member>.carries` is null where the genotype does not say (`./.`, `0/.`). `NOT carries` leaves out a member without a
  call, and `carries IS NOT TRUE` includes one.
- `<member>.zygosity` is `hom_ref`, `het` or `hom_alt`, and null wherever the call does not determine it: where
  `<member>.carries` is null, for `./1`, which carries, and in `sv` for a gain of copies. `NOT carries` and
  `zygosity = 'hom_ref'` match the same calls. `carries IS TRUE` and `zygosity IN ('het', 'hom_alt')` do not. A haploid
  call of the allele, as on a male's X, is `hom_alt`.
- In a `cnv` row, `sv.<member>.gt` is `./.` for a member at the expected copy number and `./1` for any gain. Carrying
  and zygosity come from `sv.<member>.cn` there.
- `variants.<member>.pl` keeps the likelihoods of the reference and this allele after the split, not normalised again:
  the smallest is not 0 where the member's call was of another allele of the site.
- `repeats.<member>.gt` numbers the alleles of the member's own file. Compare `repeats.<member>.repcn` between members.
- A field that the caller did not write is null, a list too: `<member>.ft` is `['PASS']` for a call that passed and null
  without FT. Only `variants.filters` and `sv.filters` can be empty lists, for a FILTER of `.`.

## Consequences and genes

- `variants.worst_consequence` and `variants.worst_impact` are the allele's own: the most severe consequence over the
  transcripts of all its genes, and that transcript's impact.
- `variants.genes` has one entry per gene, from the most to the least severe consequence, and the first is the gene of
  `variants.worst_consequence`. Group by gene through `unnest(genes)`, `variants.hgnc_ids` or `transcripts`: `genes[1]`
  alone files an allele under one of its genes.
- `variants.genes.ensembl_gene_id` is Ensembl's id, also where the entry's transcript is RefSeq's. `transcripts.gene_id`
  is the id of the transcript's own gene model, NCBI's for RefSeq. The two do not join.
- `variants.hgnc_ids` and `sv.hgnc_ids` are null where no gene of the row has an HGNC id, about one row in six. Such a
  gene is found by `transcripts.gene_symbol`.
- In `sv_transcripts.consequences` every transcript that a record touches has `feature_truncation` or
  `feature_elongation`, also within an intron. What is hit is in `sv_transcripts.complete_overlap`,
  `sv_transcripts.protein_pos`, `sv_transcripts.exon` and `sv_transcripts.intron`.

## Words and scales of other tools

Passed through unchanged: from the Annotator `variants.variant_type`, `transcripts.biotype`, `transcripts.impact`, the
consequence terms, the HGVS names (`p.(Ser134Leu)`, with parentheses) and ClinVar's texts; from the caller
`variants.filters`, `variants.qual`, `<member>.ft`, `<member>.dn`, `<member>.gq` and `<member>.pl`. `SELECT DISTINCT`
lists a column's words.

- Impact is lowercase: `high`, `moderate`, `low`, `modifier`.
- `transcripts.biotype` and `sv_transcripts.biotype` are the Annotator's words, not Ensembl's: `mRNA` where VEP has
  `protein_coding`, and `lncRNA`, `pseudogenic_transcript`, `primary_transcript`.
- `sv.sv_type` is the caller's word (`DEL`, `INS`, `DUP`, `INV`, `BND`, `CNV`), and `sv_transcripts.variant_type` the
  Annotator's for the same allele (`deletion`, `copy_number_loss`, `tandem_duplication`). A `CNV` is a segment that one
  member loses and another gains: it has rows of both types in `sv_transcripts`, and its summary columns in `sv` are
  null.
- `variants.clinvar_classification` is ClinVar's text, such as `Pathogenic/Likely pathogenic` or
  `Pathogenic; risk factor`, so match it by a pattern. It is the best-reviewed record's, and records that disagree are
  only in `variants.clinvar_records`.
- DRAGEN's GQ includes a prior against an allele that the population does not support, so a GQ floor tuned on GATK's
  calls cuts true de novo calls.

## Frequencies and intervals

- An allele that a source lacks has a null frequency, which `af < x` drops: use `coalesce`.
- The `variants.gnomad_joint_*` columns are gnomAD's own joint values, matched by position, REF and ALT exactly. Null
  means that the archive lacks the allele, which holds for every `chrM` allele and for an indel that gnomAD writes at
  another position. 0 means that gnomAD lists the allele with a count of 0, and `variants.gnomad_joint_faf95_max` is
  also 0 where gnomAD computed none.
- `sv.gnomad_sv_af` is of the best match: the gnomAD entry of the same direction with the highest reciprocal overlap,
  from 0.5 up. `sv.gnomad_sv_max_af` is the highest over all matches, and the one for a filter on common variants. A
  record with `sv.annotated_as = 'small_variant'` has neither, and has `sv.gnomad_genomes_af` and
  `sv.gnomad_exomes_af` instead.
- For any other record, an entry of `sv.clinvar_records` is of another variant that overlaps it, by
  `sv.clinvar_records.reciprocal_overlap`.
- `sv.pos` is the base before a record's bases, which run to `sv.end_pos`.
