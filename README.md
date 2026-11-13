# case-annotation-pipeline

Annotates the variant calls of one case into Parquet tables that an analyst or an agent can query directly. A case is one family: its members and everything sequenced for them.

The pipeline takes a sheet with one row for each family. A row names the family's joint VCF and its pedigree file and, where the family has them, its VCFs of structural variants, copy number segments and tandem repeats. For each family the pipeline writes one folder of tables, with a manifest that describes them.

The workflows, the tools and their documentation follow in the next changes.
