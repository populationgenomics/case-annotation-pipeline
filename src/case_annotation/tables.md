# Family tables `0.5.0`: what the tables do not show

These notes hold what a query on the tables gets wrong without an error. A name is written with its table, as in
`variants.genes.ensembl_gene_id`, and `<member>` stands for a member's column, which is named as the sample. Only the
manifest says who the proband is, and it lists the columns that the run left out.

## Rows

- `variants` holds the family's joint calls and is the place to search, with a row only for an allele that a member
  carries. A member's genome VCF adds about one call in 8,000 to it, nearly all with GQ below 20 or in conflict with
  inheritance, so reading it through or comparing it with the tables finds nothing that a search needs. Open it at one
  position only, to ask why a specific allele that was expected has no row.
- `variants` has one row per alternate allele, after `bcftools norm` split a site of several alleles and left-aligned
  its indels. `variants.pos`, `variants.ref`, `variants.alt` and `variants.variant_id` are of that form, so another
  file's spelling of an indel finds nothing.
- `transcripts.is_mane_select` and `sv_transcripts.is_mane_select` are true on RefSeq rows only, also for a transcript
  that Ensembl has.
- `sv` holds the records of two callers. Both can report one event, and nothing links the two rows. A record that is
  also a row of `variants` names it in `sv.small_variant_id`, so a query over both tables keeps
  `sv.small_variant_id IS NULL`.
- `repeats` has a row for every repeat of the catalogue, whatever the members' alleles: a row is not a finding.

## A member's call

- `variants.<member>.gt` is the genotype after the split: `1/2` is `1/0` in one row and `0/1` in the other, and `het` in
  both. A `0`, `hom_ref` and a false `<member>.carries` say "not this allele": another allele of the site reads as `0`.
- `variants.<member>.allele_balance` compares the allele's reads with the reference's alone, so it is near 1 in both
  rows of a `1/2` call.
- `<member>.carries` is null where the genotype does not say (`./.`, `0/.`), so a false one is a member who was called.
- `<member>.zygosity` is null wherever the call does not determine it: where `<member>.carries` is null, for `./1`,
  which carries, and in `sv` for a gain of copies. A haploid call of the allele is `hom_alt`.
- In a row of the copy-number caller, `sv.<member>.gt` does not say who carries: it can read `./.` at the expected copy
  number, and names one allele for a gain of any size. `sv.<member>.cn` decides carrying and zygosity there.
- A field that one caller does not write is null in all its rows of `sv`, as `sv.<member>.gq` can be in a copy-number
  row, whose quality is `sv.<member>.qs`: a floor on it removes every such row.
- `repeats.<member>.gt` numbers the alleles of the member's own file. Compare `repeats.<member>.repcn` between members.
- `variants.filters`, `sv.filters`, `repeats.<member>.filters` and `<member>.ft` hold `['PASS']` for a pass. The first
  three are never null, and empty for a FILTER of `.`. `<member>.ft` is never empty, and null for a call without FT.
- `<member>.dn` and `<member>.dq` are null in every row of a member that the caller did not call as a child with its
  parents: a comparison on them matches no row.
- Some callers fold a population prior into `<member>.gq`, which lowers it for a true call that the population does not
  support, as a newly arisen allele is, so a floor on it can drop those calls. Such a caller's `<member>.pl` can hold
  the likelihoods without that prior.
- The caller's own words and scales are passed through, and no column or manifest key names the caller: the names in
  those filter lists, the values of `<member>.dn` and `repeats.<member>.so`, and the scales of `variants.qual`,
  `sv.qual`, `<member>.gq`, `<member>.pl`, `<member>.dq`, `variants.<member>.msc`, `sv.<member>.qs` and
  `sv.<member>.sm`. A literal or a threshold from another caller's habits is wrong without an error.

## Genes and consequences

- `variants.worst_consequence` and `variants.worst_impact` are over the transcripts of all the allele's genes.
  `variants.genes` has one entry per gene, the most severe first, so `variants.genes[1]` alone files an allele under one
  of its genes. Group by gene through `unnest(variants.genes)`, `variants.hgnc_ids` or `transcripts`.
- An allele whose only hit is an exon of a non-coding RNA gene has the impact `modifier`, so no filter on
  `variants.worst_impact` finds it. Its consequence does: `non_coding_transcript_exon_variant` in
  `variants.genes.worst_consequence`.
- `variants.genes.ensembl_gene_id` is Ensembl's id, also where the entry's transcript is RefSeq's. `transcripts.gene_id`
  is NCBI's on a RefSeq row, so a join of the two loses the RefSeq rows.
- `variants.hgnc_ids` and `sv.hgnc_ids` are null for a row without a transcript, and also where no gene of the row has
  an HGNC id. Such a gene is found by `transcripts.gene_symbol` or `sv_transcripts.gene_symbol`.
- In `sv_transcripts.consequences`, every transcript that a deletion, an insertion or a duplication touches has
  `feature_truncation` or `feature_elongation`, also within an intron. What is hit is in
  `sv_transcripts.complete_overlap`, `sv_transcripts.protein_pos`, `sv_transcripts.exon` and `sv_transcripts.intron`.
- The words are the Annotator's, and a wrong one matches nothing. `transcripts.impact` is lowercase: `high`, `moderate`,
  `low`, `modifier`. `transcripts.biotype` has `mRNA` where VEP has `protein_coding`. The consequence terms are VEP's
  and some of the Annotator's own, such as `inframe_indel`. A protein's HGVS name has parentheses: `p.(Ser134Leu)`.
- `sv.sv_type` has `DEL` where `sv_transcripts.variant_type` has the Annotator's `deletion` or `copy_number_loss`. A
  `CNV` there is a segment that one member loses and another gains: it has rows of both types in `sv_transcripts`, its
  best-match and frequency columns in `sv` are null, its matches are in `sv.gnomad_sv_matches`, and
  `sv.<member>.carries` is true for a loss and for a gain.

## Frequencies and ClinVar

- An allele that a source lacks has a null frequency, which a comparison drops.
- The `variants.gnomad_joint_*` columns are matched to gnomAD's joint release by position, REF and ALT exactly. Null
  means that the archive lacks the allele, which holds for every `chrM` allele and for an indel that gnomAD writes at
  another position. 0 means a count of 0 in gnomAD, and `variants.gnomad_joint_faf95_max` is also 0 where gnomAD
  computed none.
- For a deletion, a duplication, an inversion or a copy number segment, `sv.gnomad_sv_af` is of the gnomAD entry of the
  same direction with the highest reciprocal overlap, from 0.5 up. `sv.gnomad_sv_max_af` is the highest over all such
  entries, those that gnomAD filtered included, and is the one for a filter on common variants.
- For an insertion, the entries are those within about 500 bases of `sv.pos`, of any length, and `sv.gnomad_sv_af` is of
  the nearest (`sv.gnomad_sv_distance`): an insertion is known nearby, which does not say that this one is common.
- A record with `sv.annotated_as = 'small_variant'` has neither. It has `sv.gnomad_genomes_af`, `sv.gnomad_exomes_af`
  and the ClinVar records of its own allele instead.
- A breakend and an insertion not annotated as a small variant are compared with no ClinVar record: a null
  `sv.clinvar_classification` there does not say that ClinVar lacks them.
- `variants.clinvar_classification` is the best-reviewed record's, and records that disagree are only in
  `variants.clinvar_records`. For any other record of `sv`, `sv.clinvar_classification` is the best-overlapping
  record's, of another variant that overlaps it at 0.5 or more (`sv.clinvar_overlap`). A null there does not say that
  nothing pathogenic lies inside a large deletion.
- A classification is ClinVar's text, such as `Pathogenic/Likely pathogenic`. A pattern for it has to leave out
  `Conflicting classifications of pathogenicity`.
