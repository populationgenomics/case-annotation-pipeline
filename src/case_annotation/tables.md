# Family tables `0.5.1`: what the tables do not show

Only the manifest says who the proband is.

## Rows

- `variants` has one row per alternate allele that a member carries, after `bcftools norm` (split, left-aligned). In the
  one trio compared, a member's genome VCF held about one call in 8,000 more, nine in ten of them with GQ below 20.
- `transcripts.is_mane_select` and `sv_transcripts.is_mane_select` are true on RefSeq rows only.
- `sv` holds two callers' records: one event can be a row of each, and nothing links the two.
- `repeats` has a row for every repeat of the catalogue, whatever the alleles.
- `roh.chrom` and `upd.chrom` are chr1 to chr22 only.
- `upd` has small regions in every trio (84 in the one ordinary trio measured, 80 of them under 100 kb), and no row for
  a child in the manifest's `upd_not_assessed`.

## A member's call

- A `1/2` call is `1/0` in one row of `variants` and `0/1` in the other, and `het` in both, so a `0` in
  `variants.<member>.gt` or a `hom_ref` can be another allele of the site.
- `variants.<member>.allele_balance` is the allele's reads over the allele's plus the reference's, so near 1 in both
  rows of a `1/2` call.
- `<member>.carries` is null for `./.` and `0/.`. `<member>.zygosity` is also null for `./1` and for a gain of copies in
  `sv`, and `hom_alt` for a haploid call of the allele.
- In a copy-number row, `sv.<member>.gt` can be `./.` at the expected copy number and has one allele for a gain of any
  size. `<member>.carries` and `<member>.zygosity` come from `sv.<member>.cn` there.
- `repeats.<member>.gt` numbers the alleles of the member's own file and does not compare between members.
- `<member>.dn` and `<member>.dq` are null for a member that the caller did not call as a child with its parents.
- Some callers fold a population prior into `<member>.gq`, which lowers it for a true call that the population does not
  support, such as a newly arisen allele. Their `<member>.pl` can be without the prior.

## Genes and consequences

- `transcripts.impact`, `variants.worst_impact` (over all the allele's genes) and the consequence terms are the
  Annotator's and take no account of the columns that other steps add: `transcripts.utr5_*`, `transcripts.nmd_*`,
  `transcripts.ag_gains`, `variants.branch_points`, `variants.gencode_polya_signal`, `variants.aparent2_*` and
  AlphaGenome's `variants.avi_*`.
- An allele in an exon of a non-coding RNA gene can have the impact `modifier` alone
  (`non_coding_transcript_exon_variant`).
- The consequence terms include some that VEP lacks, such as `inframe_indel`.
- `variants.genes.ensembl_gene_id` is Ensembl's, also where the entry's transcript is RefSeq's. `transcripts.gene_id` is
  NCBI's on a RefSeq row.
- `variants.hgnc_ids` and `sv.hgnc_ids` are null also where a row has genes and none has an HGNC id.
- `sv_transcripts.consequences` has `feature_truncation` or `feature_elongation` for every transcript that a deletion,
  an insertion or a duplication touches, also within an intron.
- A `CNV` in `sv.sv_type` is a segment that one member loses and another gains: its best-match and frequency columns are
  null, its matches are in `sv.gnomad_sv_matches` and `sv.clinvar_records`, and `sv.<member>.carries` is true for a loss
  and for a gain.

## Frequencies and ClinVar

- The `variants.gnomad_joint_*` columns are matched by position, REF and ALT exactly, so they are null for an indel that
  gnomAD writes at another position.
- `sv.gnomad_sv_af` is of the gnomAD entry of the same direction with the highest reciprocal overlap, from 0.5 up. For
  an insertion it is of the nearest entry within about 500 bases of `sv.pos`, of any length. `sv.gnomad_sv_max_af` is
  the highest over all matched entries, those that gnomAD filtered included.
- A record with `sv.annotated_as = 'small_variant'` has neither. It has `sv.gnomad_genomes_af`, `sv.gnomad_exomes_af`
  and the ClinVar records of its own allele.
- A breakend and an insertion not annotated as a small variant never hold a ClinVar record.
- `variants.clinvar_classification` is the best-reviewed record's, and records that disagree are only in
  `variants.clinvar_records`. For any other record of `sv`, `sv.clinvar_classification` is that of the best-overlapping
  ClinVar variant, at an overlap of 0.5 or more.
