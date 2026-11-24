# Glossary

This page defines the words that the other documents use, each in a sentence or two, with a link to the place that
explains it. The documents use one word for one thing, so a word means here what it means everywhere else in the
repository.

The entries are grouped by topic. The examples use the same invented data as every other document: the family `family-a`
with the samples `S0101` to `S0103`. Every id and value in the examples is invented.

## Genetics and sequencing

A reader who works with sequencing data can skip this group.

### Family

The unit that the pipeline works on: the people of one pedigree, with the files that hold what was sequenced for them. A
family has an id, such as `family-a`, which names its row of the families sheet, the family column of its PED file and
its folder of tables. A case, in the sense of one analysis, is one family.

### Proband

The person whose disease the analysis of a family is about, usually an affected child. The other family members are
described by their relation to the proband. A families sheet names each family's proband, or none
([The pedigree and the roles](family-tables.md#the-pedigree-and-the-roles)).

### Trio

A proband and both parents, all three sequenced.

### Pedigree

Who is whose parent in a family, with each person's sex and whether they are affected. A PED file writes it as one line
for each person with six columns: family, individual, father, mother, sex and affected status.

### Reference genome

The agreed sequence that reads are compared against. Everything in this repository uses GRCh38, with chromosomes named
`chr1` to `chr22`, `chrX`, `chrY` and `chrM`.

### Variant

A place where a person's DNA differs from the reference genome. A single-base change is an SNV, a short insertion or
deletion an indel, and larger changes are structural variants and copy number changes. An allele is one version of the
sequence at a position: the reference allele, or an alternate one.

### Chromosome copy

A person has two copies of each chromosome, one from each parent, so a variant sits on one copy or on both. Males have
one copy of chrX outside the pseudoautosomal regions, and everyone has one copy of chrM.

### Genotype

Which alleles a person carries at a position, written as in a VCF: `0/0` for two reference alleles, `0/1` for one
alternate allele (heterozygous), `1/1` for two (homozygous), and a single `1` where the person has one copy
(hemizygous).

### De novo variant

A variant that a child carries and neither parent does, because it arose new in the child
([De novo variants](annotation.md#de-novo-variants)).

### Upstream reading frame

A stretch of a transcript's 5′ untranslated region (5′ UTR), the part before the protein's own start codon, that begins
with a start codon of its own. A ribosome that starts there can fail to reach the protein's start, so a variant that
creates or removes such a frame can change how much protein is made. A frame that ends before the protein's start is a
uORF, and one that runs past it is an overlapping one, an oORF
([Upstream reading frames and nonsense-mediated decay](annotation.md#upstream-reading-frames-and-nonsense-mediated-decay)).

### Nonsense-mediated decay

How a cell destroys a transcript whose reading frame ends early, so that a variant which stops a protein short usually
leaves no protein at all. A variant near the transcript's end, or near its start, can escape the decay, and the
shortened protein is then made.

### Branch point

The base inside an intron, most often an adenine 18 to 44 bases before the intron's end, that the cell joins the start
of the intron to when it splices the intron out. A variant that changes it, or the short motif around it, can leave the
intron in the transcript or drop the next exon ([Catalogued branch points](family-tables.md#catalogued-branch-points)).
The cell takes the first `AG` behind the branch point as the intron's end, so a variant that creates an `AG` between the
two can move that end ([AG gained before an acceptor](family-tables.md#ag-gained-before-an-acceptor)).

### Structural variant

A change of 50 bases or more, or a rearrangement: a deletion, an insertion, a duplication, an inversion, or a junction
between two distant places. A VCF writes a junction as two breakends, one record for each side, which name each other as
mates. DRAGEN's structural-variant caller finds them from reads that span or are split by a junction, and also reports
some shorter insertions and deletions, which the small-variant caller reports too. In the family tables its records are
rows of `sv`, beside the copy-number caller's
([Structural variants, copy number and repeats](family-tables.md#structural-variants-copy-number-and-repeats)).

### Copy number segment

A stretch of a chromosome over which a person has fewer or more copies than expected, which is two, or one on a male's X
and Y outside the pseudoautosomal regions. A copy-number caller finds segments from the depth of the reads and gives
each member's copy number of each. The family tables hold the copy-number variants in `sv`, beside the structural
variants and told apart by the `caller` column: a row is a segment whose copy number differs from the expected one, and
the segments at the expected copy number are left out. A member carries it by copy number, a gain with no zygosity, and
a deletion that both callers find has a row from each.

### Tandem repeat

A short unit of sequence, such as `CAG`, repeated end to end. Some diseases come from a repeat that has grown beyond a
size, a repeat expansion. A repeat genotyper sizes a catalogue of known repeats in one sample, in repeat units for each
of the sample's alleles.

### Run of homozygosity

A stretch of a chromosome in which a person's two copies are the same, so that every site in it is homozygous. A long
one arises where both copies come from one ancestor, as in the child of related parents, where one copy is deleted, and
in an isodisomy ([Uniparental disomy](#uniparental-disomy)). The table `roh` holds the runs that `bcftools roh` reports
([Runs of homozygosity and uniparental disomy](family-tables.md#runs-of-homozygosity-and-uniparental-disomy)).

### Uniparental disomy

Both copies of a chromosome, or of a part of one, from the same parent, and none from the other. In a heterodisomy the
child has both of that parent's copies, and in an isodisomy one of them twice, which is then a run of homozygosity. A
deletion of one parent's copy looks the same in the genotypes as an isodisomy of the other parent's. The table `upd`
holds the regions that the tool `upd` reports for a child whose parents are both sequenced.

### Pseudoautosomal regions

The two stretches at the ends of chrX and chrY that are alike on both, written PARs. Inside them every person has two
copies, as on an autosome.

### Variant calling

Comparing a person's reads with the reference genome and writing down the variants. A VCF file lists variants with each
sample's genotype, a gVCF also records the stretches where a sample matches the reference, and a joint VCF holds the
calls made over a whole family at once.

### Sample name

The name that a family's VCF gives one member's sample, such as `S0103`. It is the member's individual in the PED file,
the member's column in `variants` and the member's `sample_id` in the manifest.

### DRAGEN

Illumina's software that maps reads and calls variants. The pipeline runs none of its callers. It reads VCFs in the
layouts that DRAGEN's callers write ([The input contract](annotation.md#the-input-contract)), and it runs the Annotator,
which ships inside DRAGEN Software Mode.

### Personalisation

A DRAGEN step that infers a sample's own haplotypes from the whole genome's reads and uses them when it calls. With it
on, the genotype likelihoods (PL) of a known and a novel allele are not comparable
([De novo variants](annotation.md#de-novo-variants)).

### ClinVar and gnomAD

Two public databases that an analysis leans on. ClinVar collects laboratories' classifications of variants, and gnomAD
gives how often each variant occurs in the population.

## Annotation

### Annotation

The facts that an analysis adds to each variant before anyone reads it: its effect on each transcript, its frequency in
the population, prediction scores and its ClinVar classification ([Annotation](annotation.md)).

### Annotator

Illumina's standalone DRAGEN Annotator, which `annotator.nf` runs. Its data and its output are licensed to the
organisation that runs it, and the pipeline takes its image as a param ([The Annotator's image](annotator-image.md)).

### Gene models

Where each gene's transcripts lie on the genome: their exons and their coding sequence. Ensembl and RefSeq each publish
a set, in releases, and a transcript's id carries a version that changes with its model. The Annotator reads both sets
from its catalog, and the rules of VEP's two plugins read the same releases from two prepared files
([The gene models](annotation.md#the-gene-models)).

### Catalog

Illumina's word, in Illumina's spelling, for a dated set of data source versions, with an id such as `202610`. It has
nothing to do with the catalogue of repeats that a repeat genotyper sizes.

### Dated source

An annotation source that carries clinical assertions or literature, such as ClinVar or OMIM, so that what it says
depends on the date of its version. Every other source is a current source. A run includes a dated source only by naming
it, and a family's manifest records the versions that the run pinned and the sources that it left out
([Dated sources](annotation.md#dated-sources)). STRchive is a dated source too, which a run includes by giving the file
of one of its releases.

### Blinded positions

Positions that a run is given for one family, at which the family's tables hold a null in the ClinVar columns, as for a
variant without a record. A row of `variants` is at such a position where its footprint, the bases of its `ref` and for
an insertion the base after them, shares a base with one of the positions that were given for a small variant. The
unions keep the values. Nothing in the family's folder says that the run was given such positions, or where they are
([Blinded positions](annotation.md#blinded-positions)).

### STRchive

A curated, public catalogue of the tandem repeats that are linked to a disease, with the sizes at which each repeat
counts as benign, intermediate and pathogenic. A run that is given the loci of one of its releases adds those sizes to
`repeats` ([STRchive's ranges of the repeats](annotation.md#strchives-ranges-of-the-repeats)).

### AlphaGenome

Google DeepMind's model of DNA sequence. DeepMind publishes tables of scores that it computed with the model for every
single-base change, among them the AlphaGenome Variant Impact (AVI) score. An echtvar archive can be built from each
table ([AlphaGenome's archives](annotation.md#alphagenomes-archives)), and a run that is given one adds the scores to
`variants` ([AlphaGenome's scores](family-tables.md#alphagenomes-scores)).

### Families sheet

The input of one run of `family_tables.nf`: a CSV with one row for each family of a group, naming its joint VCF, its PED
file and its proband, and the VCFs of its larger variants where it has them
([The input contract](annotation.md#the-input-contract)).

### Union and group

A group is the set of families that one run of `family_tables.nf` annotates together, and the union is the set of all
their alleles, which is annotated once ([Family tables](annotation.md#family-tables)). The structural variants and copy
number segments of a group have a union of their own, annotated in a second run.

### Family tables

The two Parquet tables and the manifest that `family_tables.nf` writes for one family: `variants`, which holds each
allele with every sequenced member's call of it, and `transcripts`. They answer by region, gene or variant in well under
a second ([Reading the family tables](family-tables.md)). A family that also brings the calls of its larger variants has
up to three more tables. `sv` holds the structural variants and the copy-number variants, the records of two callers,
told apart by its `caller` column; `sv_transcripts` holds the transcripts that those records hit; and `repeats` holds
the tandem repeats. `sv` and `repeats` have a column for every sequenced member as well. A copy-number row of `sv` is a
segment whose copy number differs from the expected one. A member carries it by copy number, a gain with no zygosity,
and one event can have a row from each caller. Every family also has `roh`, with each sequenced member's runs of
homozygosity, and a family with a child whose parents are both sequenced has `upd`, with that child's regions of
uniparental disomy.

### Manifest

The JSON index of a family's folder, `manifest.json`: the tables that the folder holds with their sizes and checksums,
the pedigree with each member's role, and what annotated the tables ([The manifest](family-tables.md#the-manifest)).

### Larger variants

What the documents call the variants that a joint VCF does not hold and that other callers write into VCFs of their own:
structural variants, copy number segments and tandem repeats. The family tables hold them as the callers wrote them, the
structural variants and copy number segments with the Annotator's annotation and the repeats with STRchive's size ranges
where the run was given a release of STRchive
([Structural variants, copy number and repeats](family-tables.md#structural-variants-copy-number-and-repeats)).

### Reciprocal overlap

How alike two intervals of a chromosome are: the bases they share, divided by the length of the longer one. It is 1 for
the same interval, and 0.5 where an interval of 6 kb lies inside one of 12 kb. The Annotator matches a structural
variant to a database's entry at a reciprocal overlap of 0.5 or more
([What the annotation says of a record](family-tables.md#what-the-annotation-says-of-a-record)).

### Variant id

The key that joins the family tables, of the form `<chrom>-<pos>-<ref>-<alt>`, such as `chr2-1234567-C-T`. It names one
alternate allele.
