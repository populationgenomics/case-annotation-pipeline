"""An invented file of STRchive's loci, in the layout of the file that a release of STRchive publishes.

Every locus, gene, range and release here is invented. The loci are placed on the invented catalogue of repeats of
structural_vcfs, whose repeats lie at chr2:700-721 (GCC), chr1:2000-2030 (CAG) and chrX:9000-9060 (CGG).
"""

from __future__ import annotations

import datetime as dt
import json
import pathlib
from collections.abc import Sequence
from typing import Any

from case_annotation import strchive

# The release that a run states for the file, and the day on which it says the release was published.
RELEASE = 'v0.0.0'
RELEASE_DATE = dt.date(2099, 1, 1)
FILE_NAME = 'STRchive-loci-v0.0.0.invented.json'


def locus(**changed: Any) -> dict[str, Any]:  # noqa: ANN401 (JSON values)
    """An entry of the file with the fields that the tables read, and one that they do not; `changed` replaces some."""
    entry: dict[str, Any] = {
        'id': 'INV1_GENEA',
        'gene': 'GENEA',
        'evidence': ['Definitive'],
        'chrom': 'chr1',
        'start_hg38': 1_990,
        'stop_hg38': 2_040,
        'inheritance': ['AD'],
        'reference_motif_reference_orientation': ['AGC'],
        'pathogenic_motif_reference_orientation': ['AGC'],
        'benign_min': 5,
        'benign_max': 20,
        'intermediate_min': 21,
        'intermediate_max': 35,
        'pathogenic_min': 36,
        'pathogenic_max': 200,
        'details': 'Invented.',
    }
    return entry | changed


# One locus of the repeat on chr1, with a wider interval and another rotation of its unit. One of the repeat on chr2,
# whose motifs hold the code N, whose benign range is one size and whose link to its disease is disputed. One that
# begins at the base after the repeat on chrX, and so is of no repeat, and that has no level of evidence.
LOCI: Sequence[dict[str, Any]] = (
    locus(),
    locus(
        id='INV2_GENEB',
        gene='GENEB',
        evidence=['Disputed'],
        chrom='chr2',
        start_hg38=695,
        stop_hg38=725,
        inheritance=['AD', 'AR'],
        reference_motif_reference_orientation=['GCN'],
        pathogenic_motif_reference_orientation=['NGC'],
        benign_min=7,
        benign_max=7,
        intermediate_min=None,
        intermediate_max=None,
        pathogenic_min=9,
        pathogenic_max=14,
    ),
    locus(
        id='INVX_GENEX',
        gene='GENEX',
        evidence=[],
        chrom='chrX',
        start_hg38=9_060,
        stop_hg38=9_100,
        inheritance=[],
        reference_motif_reference_orientation=['CGG'],
        pathogenic_motif_reference_orientation=['CGG'],
    ),
)
# The repeats of structural_vcfs that have a locus here, each with its locus, and the one that has none.
MATCHED = {'RPTA': 'INV1_GENEA', 'RPTB_GCC': 'INV2_GENEB'}
UNMATCHED = ('RPTX',)


def write(path: pathlib.Path, loci: object = LOCI) -> pathlib.Path:
    """Writes the loci, or anything else that JSON can hold, as a release's file."""
    path.write_text(json.dumps(loci, indent=2) + '\n')
    return path


def release(folder: pathlib.Path, loci: object = LOCI) -> strchive.Release:
    """The invented release as a run gives it: its file in the folder, with its name and date as stated."""
    return strchive.Release(file=write(folder / FILE_NAME, loci), name=RELEASE, date=RELEASE_DATE)
