"""
Alignment statistics for Stockholm files.

Quality metrics for evaluating multiple sequence alignments.
"""


def compute_pairwise_identity(sequence_entries):
    """
    Compute average pairwise identity for each sequence in the alignment.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples (aligned, with gaps)

    Returns:
        list: [(seq_name, avg_identity, num_comparisons)] sorted by avg_identity ascending
    """
    gaps = set('.-')
    n = len(sequence_entries)
    if n < 2:
        return []

    results = []
    for i, (name_i, seq_i) in enumerate(sequence_entries):
        identities = []
        for j, (name_j, seq_j) in enumerate(sequence_entries):
            if i == j:
                continue
            matches = 0
            compared = 0
            for ci, cj in zip(seq_i, seq_j):
                if ci in gaps or cj in gaps:
                    continue
                compared += 1
                if ci.upper() == cj.upper():
                    matches += 1
            if compared > 0:
                identities.append(matches / compared * 100)
        avg = sum(identities) / len(identities) if identities else 0
        results.append((name_i, avg, len(identities)))

    results.sort(key=lambda x: x[1])
    return results
