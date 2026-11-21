# Notices of the works that this package restates

Three modules of this package restate, in Python, rules that others published as Perl. No file of theirs is copied. The
licences of those works ask that their notices stay with what is made from them, so they are kept here and are packed
with the package.

| Module                              | Restates                                                                                                                                                                                                                                                                                                | Licence    | Notice                                               |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- | ---------------------------------------------------- |
| `utr5`                              | the UTRannotator plugin for Ensembl VEP, [ImperialCardioGenetics/UTRannotator](https://github.com/ImperialCardioGenetics/UTRannotator) at the commit `3d9d466e356752aa69d61b62da493edddf374363`: `UTRannotator.pm` and the six files under `src/`                                                       | MIT        | [UTRannotator-LICENSE.txt](UTRannotator-LICENSE.txt) |
| `nmd`                               | the NMD plugin, `NMD.pm` of [Ensembl/VEP_plugins](https://github.com/Ensembl/VEP_plugins) on the branch `release/116` at the commit `2341cdc2c173caeae043e60de412fa185efe00bc`                                                                                                                          | Apache 2.0 | [ENSEMBL-NOTICE.txt](ENSEMBL-NOTICE.txt)             |
| `plugin_gates`, `transcript_models` | how Ensembl VEP 116.3 reads a GFF file into transcripts and gives an allele its consequence terms: modules of [ensembl-vep](https://github.com/Ensembl/ensembl-vep), [ensembl-variation](https://github.com/Ensembl/ensembl-variation) and [ensembl](https://github.com/Ensembl/ensembl) at release 116 | Apache 2.0 | [ENSEMBL-NOTICE.txt](ENSEMBL-NOTICE.txt)             |

The Apache License 2.0 is this repository's own licence, and its text is the file `LICENSE` at the repository's root.
Each of the modules says in its opening lines what it restates and where it departs from it, which is the statement of
changes that the licence asks for.

[vep-rs](https://github.com/natera-open-source/vep-rs) (Apache 2.0, Copyright (c) 2026-present Natera, Inc.) was read at
its tag `v0.6.0` for the logic of both plugins and for its account of the defects of Ensembl's edition of UTRannotator.
Nothing of it is copied or restated.
