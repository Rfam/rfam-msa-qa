"""
Configuration parameters for Stockholm validation and fixing.

Modify these values to adjust the behavior of the validation script.
"""

# =============================================================================
# BLAST Search Parameters
# =============================================================================

# Minimum identity percentage to accept a BLAST hit (0-100)
BLAST_MIN_IDENTITY = 90 #95

# Minimum query coverage percentage to accept a BLAST hit (0-100)
BLAST_MIN_COVERAGE = 80 #90

# Maximum e-value to accept a BLAST hit
BLAST_MAX_EVALUE = 1e-5 #1e-10

# Number of BLAST hits to retrieve (top N hits)
BLAST_HITLIST_SIZE = 5


# =============================================================================
# NCBI Request Parameters
# =============================================================================

# Delay between NCBI requests in seconds (to avoid rate limiting)
NCBI_REQUEST_DELAY = 0.5
