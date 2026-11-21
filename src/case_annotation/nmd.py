"""Whether an allele that ends a reading frame early is predicted to escape nonsense-mediated decay.

The columns are described in docs/annotation.md.

This is the rule of the NMD plugin of Ensembl's collection for VEP (https://github.com/Ensembl/VEP_plugins, Apache
2.0, licenses/ENSEMBL-NOTICE.txt), restated in Python from `NMD.pm` on the branch `release/116` at the commit
SOURCE_COMMIT. The rule of vep-rs (https://github.com/natera-open-source/vep-rs, Apache 2.0, tag v0.6.0,
`crates/vep-builtins/src/plugins/nmd.rs`) was read beside it. No code of either is copied.

The plugin looks at one base of the allele, its last on the genome, which for an insertion is the base before it,
and at the coding position of the allele's last base on the transcript. It predicts an escape for four reasons and
gives the first that applies, in this order:

- `intronless`: the transcript has one exon;
- `E-1`: the base lies in the transcript's last exon, its untranslated part included;
- `E-2_end-50`: the base lies in the last 52 bases of the last exon but one, which the plugin calls its last 50;
- `cds_start+100`: the allele ends in the first 101 coding bases, which the plugin calls the first 100.

The numbers are the plugin's own and are kept as they are. Which alleles the plugin evaluates at all is VEP's
decision (plugin_gates).
"""

from __future__ import annotations

from case_annotation import transcript_models

# The commit of https://github.com/Ensembl/VEP_plugins that last changed NMD.pm on the branch release/116.
SOURCE_COMMIT = '2341cdc2c173caeae043e60de412fa185efe00bc'

REASONS = ('intronless', 'E-1', 'E-2_end-50', 'cds_start+100')
_PENULTIMATE_REACH = 51
_CODING_REACH = 101


def escape_reason(transcript: transcript_models.Transcript, last_base: int, coding_end: int) -> str | None:
    """Why the plugin predicts that an allele escapes nonsense-mediated decay on a transcript, or None.

    `last_base` is the allele's last base on the genome, and `coding_end` the coding position of the allele's last
    base on the transcript, counted from 1.
    """
    exons = transcript.exons
    if len(exons) == 1:
        return 'intronless'
    if exons[-1][0] <= last_base <= exons[-1][1]:
        return 'E-1'
    start, end = exons[-2]
    reach = (end - _PENULTIMATE_REACH, end) if transcript.strand == 1 else (start, start + _PENULTIMATE_REACH)
    if reach[0] <= last_base <= reach[1]:
        return 'E-2_end-50'
    if coding_end <= _CODING_REACH:
        return 'cds_start+100'
    return None
