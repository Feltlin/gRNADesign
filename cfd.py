"""CFD (Cutting Frequency Determination) scoring — Doench/CRISPOR logic.

Everything else in the pipeline that needs a mismatch-tolerance number goes
through calc_cfd / mm_activity. Do not hand-roll mismatch keys elsewhere.
"""

import pickle

_BASES = ("A", "C", "G", "T")

# complement used only for building mismatch-table keys (the 'd' side is
# always DNA A/C/G/T; sg has already been through T->U, so 'U' must map
# back to 'A' here).
_KEY_COMPLEMENT = {"A": "T", "T": "A", "U": "A", "C": "G", "G": "C"}


def _complement(base):
    return _KEY_COMPLEMENT[base]


def load_cfd(mm_path, pam_path):
    """Load the mismatch and PAM CFD tables, identifying each by key shape
    (not filename) so a swapped argument order is corrected automatically."""
    with open(mm_path, "rb") as f:
        a = pickle.load(f)
    with open(pam_path, "rb") as f:
        b = pickle.load(f)

    def looks_like_mismatch(d):
        k = next(iter(d))
        return ":" in k and "," in k

    def looks_like_pam(d):
        k = next(iter(d))
        return len(k) == 2 and all(c in _BASES for c in k)

    a_is_mm, b_is_mm = looks_like_mismatch(a), looks_like_mismatch(b)
    a_is_pam, b_is_pam = looks_like_pam(a), looks_like_pam(b)

    if a_is_mm and b_is_pam:
        mm_scores, pam_scores = a, b
    elif b_is_mm and a_is_pam:
        mm_scores, pam_scores = b, a
        print("WARNING: load_cfd received mm_path/pam_path swapped; "
              "correcting based on key shape.")
    else:
        raise ValueError(
            "Could not identify mismatch/PAM tables by key shape: "
            f"first path keys look like {next(iter(a))!r}, "
            f"second path keys look like {next(iter(b))!r}"
        )

    return mm_scores, pam_scores


def calc_cfd(wt, sg, pam, mm_scores, pam_scores):
    """Canonical Doench/CRISPOR CFD score.

    wt: intended/on-target 20-mer protospacer, 5'->3', DNA letters.
    sg: compared 20-mer (variant or off-target protospacer), same frame.
    pam: 2-nt key into pam_scores (the "GG" of the NGG).
    """
    score = 1.0
    wt = wt.replace("T", "U")
    sg = sg.replace("T", "U")
    for i in range(20):
        if wt[i] == sg[i]:
            continue
        key = "r" + wt[i] + ":d" + _complement(sg[i]) + "," + str(i + 1)
        score *= mm_scores[key]
    score *= pam_scores[pam]
    return score


def mm_activity(pos0, guide_base, target_base, mm_scores, pam_scores):
    """Single-position mismatch activity (pos0 is 0-indexed spacer position).

    1.0 if bases match. A gap at target_base means an indel at that column,
    which frameshifts the guide, so activity is 0 (not looked up in the CFD
    table, which has no notion of indels).
    """
    if target_base == guide_base:
        return 1.0
    if target_base == "-" or guide_base == "-":
        return 0.0

    filler = "A"
    wt = [filler] * 20
    sg = [filler] * 20
    wt[pos0] = guide_base
    sg[pos0] = target_base
    raw = calc_cfd("".join(wt), "".join(sg), "GG", mm_scores, pam_scores)
    return raw / pam_scores["GG"]


def cfd_position_weights(mm_scores, pam_scores):
    """Data-independent per-position weight vector (length 20).

    w_i = 1 - mean over all ordered (x, y), x != y in {A,C,G,T} of
    mm_activity(i, x, y). Automatically seed-heavy if the CFD table is
    correctly indexed (higher weight = harsher penalty = seed region).
    """
    weights = []
    pairs = [(x, y) for x in _BASES for y in _BASES if x != y]
    for i in range(20):
        activities = [mm_activity(i, x, y, mm_scores, pam_scores) for x, y in pairs]
        mean_activity = sum(activities) / len(activities)
        weights.append(1.0 - mean_activity)
    return weights
