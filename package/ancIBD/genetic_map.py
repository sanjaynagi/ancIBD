"""Canonical lengths of the human autosomes along the genetic map.

Single source of truth for the genome-wide constants used both when fitting Ne
(ancIBD.ibd_stats.estimate_Ne) and when plotting karyotypes
(ancIBD.plot.plot_karyotype), so the two cannot drift apart.
Lengths are from the Eigenstrat 1240k genetic map.
"""

### Length of each autosome [in cM], keyed by chromosome number.
CH_LENGTHS_CM = {1: 286.279, 2: 268.840, 3: 223.361, 4: 214.688,
                 5: 204.089, 6: 192.040, 7: 187.221, 8: 168.003,
                 9: 166.359, 10: 181.144, 11: 158.219, 12: 174.679,
                 13: 125.706, 14: 120.203, 15: 141.860, 16: 134.038,
                 17: 128.491, 18: 117.709, 19: 107.734, 20: 108.267,
                 21: 62.786, 22: 74.110}


def chrom_length_cm(ch):
    """Return the genetic length of autosome ch [in cM]."""
    return CH_LENGTHS_CM[ch]


def chrom_length_morgan(ch):
    """Return the genetic length of autosome ch [in Morgan]."""
    return CH_LENGTHS_CM[ch] / 100.0
