"""The per candidate metric.

Conservation (entropy), PAM conservation and cutting tolerance (CFD) are kept
as separate column on purpose. They answer different question and must never
be fused into one number.
"""

import numpy as np

import cfd
from pam import motif_mask


def entropy_mean(entropy, column):
    r"""Plain mean Shannon entropy over the spacer column.

    \bar{H} = (1/L) \sum_{j \in spacer} H_j, in [0, 1]. Lower is more
    conserved across the panel.
    """
    return float(np.mean(entropy[column]))


def cfd_mean(matrix, candidate, motif, table, pam_score):
    r"""Mean CFD of the consensus guide against every aligned sequence here.

    \overline{CFD} = (1/N) \sum_k CFD(guide, variant_k) and, for one variant,

        CFD = PAM(variant PAM) * \prod_{i=1}^{20} a_i(guide_i, variant_i)

    where a_i is the position specific mismatch activity. A variant carrying a
    gap or an ambiguity code inside the protospacer, or lacking a usable PAM,
    scores 0 and still counts in the denominator, so the mean is over the
    whole panel.

    Returns nan when the spacer length puts it outside the published model.
    """
    if table is None:
        return float("nan")

    index = {base: position for position, base in enumerate(cfd.BASE)}
    product = np.ones(matrix.shape[0])
    for offset, column in enumerate(candidate.spacer_column):
        guide_base = candidate.spacer[offset]
        if guide_base not in index:
            return 0.0
        product *= table[offset, index[guide_base]][matrix[:, column]]

    return float(np.mean(product * pam_factor(matrix, candidate, motif, pam_score)))


def pam_factor(matrix, candidate, motif, pam_score):
    """Per variant PAM factor.

    SpCas9 NGG uses the published CFD PAM table on the GG dinucleotide. Any
    other motif falls back to a 1 or 0 motif match, because the CFD PAM table
    is SpCas9 specific and says nothing about another nuclease.
    """
    column = candidate.pam_column
    if motif.upper() == "NGG" and len(column) == 3:
        key = (matrix[:, column[1]].astype(np.int32) * 256
               + matrix[:, column[2]].astype(np.int32))
        return cfd.pam_lookup(pam_score)[key]

    mask = motif_mask(motif)
    factor = np.ones(matrix.shape[0], dtype=bool)
    for position, index in enumerate(column):
        factor &= mask[position][matrix[:, index]]
    return factor.astype(float)
