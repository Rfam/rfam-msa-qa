# rfam-msa-qa
Tools for ensuring the quality of multiple sequence alignments (MSAs) in Rfam.

## Stockholm File Validation

This repository includes a modular validation script for Stockholm format alignment files (`.so`, `.sto`, `.stk`).

### Usage

```bash
python3 validate_stockholm.py [-v] [--fix] [--output-mode {stdout,file}] <file1.so> [file2.so ...]
```

Options:
- `-v, --verbose`: Print detailed validation information
- `--fix`: Attempt to fix fixable errors automatically
- `--output-mode {stdout,file}`: Output mode for fixed files (default: file)
  - `file`: Create a new file with `_corrected` suffix
  - `stdout`: Print corrected content to stdout

### What is validated?

**Fatal Errors** (must be fixed manually):
- Missing `# STOCKHOLM 1.0` header
- Missing `//` terminator
- No sequences found in alignment
- **All sequences must have the same length**
- **Sequences must not contain whitespace characters**

**Fixable Errors** (can be auto-corrected with `--fix`):
- Duplicate sequences (same accession, coordinates, and sequence data)
- Missing coordinates (sequences without start/end positions)
- Overlapping sequences (sequences from the same accession that overlap by ≥1 bp)

**Warnings** (non-critical):
- Missing 2D structure consensus annotation (`#=GC SS_cons`)
- Lines exceeding 10,000 character limit

### Modular Architecture

The validation logic is split into separate modules in the `scripts/` directory:
- `fatal_errors.py`: Errors that cannot be automatically fixed
- `fixable_errors.py`: Errors that can be automatically corrected (including coordinate fixing and overlap removal)
- `stockholm_warnings.py`: Non-critical issues
- `parser.py`: Stockholm file parsing utilities

### Sequence Format

Sequences should follow the format: `ACCESSION/START-END` where:
- `ACCESSION` is the sequence identifier (e.g., from GenBank like `AF228364.1`)
- `START-END` are the coordinates indicating which portion of the original sequence is included (e.g., `1-74`)

Example: `AF228364.1/1-74`

Sequence data:
- May contain any characters except whitespace
- Gaps may be indicated by `.` or `-`

### Duplicate Detection and Removal

The script can detect and remove duplicate sequences using the `--fix` flag. Duplicates are defined as sequences that have:
1. The same accession/identifier
2. The same coordinates
3. The exact same sequence data

### Coordinate Fixing

When sequences are missing coordinates (e.g., `NZ_CP038662.1` instead of `NZ_CP038662.1/4723652-4723704`), the `--fix` flag will:

1. **NCBI Direct Lookup**: Download the FASTA sequence from NCBI and find coordinates by matching the alignment sequence
2. **BLAST Fallback**: If direct lookup fails, perform a BLAST search against NCBI's nucleotide database (accepts hits with ≥95% identity, ≥90% coverage, e-value ≤1e-10)
3. **Removal**: Sequences that cannot be resolved via either method are removed from the output

**Requirements**: BioPython (`pip install biopython`)

### Overlap Detection and Removal

Sequences from the same accession (species) that overlap by at least 1 bp are detected. When using `--fix`, overlapping sequences are removed using greedy interval scheduling to keep the maximum number of non-overlapping sequences.

### Sequence Validation

When using `--fix`, all sequences are validated against NCBI to ensure they match the source data at the given coordinates:

1. **NCBI Validation**: Fetches the source sequence from NCBI and compares it against the Stockholm sequence at the given coordinates
2. **BLAST Fallback**: If the accession is not found in NCBI or the sequence doesn't match, a BLAST search is performed to find the correct accession and coordinates
3. **Removal**: Sequences that fail both NCBI validation and BLAST are removed from the output

### Configuration

Parameters can be adjusted in `scripts/config.py`:
- `BLAST_MIN_IDENTITY`: Minimum identity % for BLAST hits (default: 95)
- `BLAST_MIN_COVERAGE`: Minimum coverage % for BLAST hits (default: 90)
- `BLAST_MAX_EVALUE`: Maximum e-value for BLAST hits (default: 1e-10)
- `NCBI_REQUEST_DELAY`: Delay between NCBI requests in seconds (default: 0.5)

### Examples

```bash
# Validate a single file
python3 validate_stockholm.py example_valid.so

# Validate multiple files with verbose output
python3 validate_stockholm.py -v file1.so file2.so file3.so

# Fix duplicate sequences and create corrected file
python3 validate_stockholm.py --fix file.so

# Fix and output to stdout
python3 validate_stockholm.py --fix --output-mode stdout file.so
```

### Continuous Integration

The repository includes a GitHub Action that automatically validates Stockholm files in pull requests. The CI check will:
- Trigger when a PR modifies `.so`, `.sto`, or `.stk` files
- Run the validation script on all changed files
- Pass only if all files are valid

### Stockholm Format Reference

Stockholm format is used for multiple sequence alignments. Basic structure:
```
# STOCKHOLM 1.0
AF228364.1/1-74    CGGCAGAUGAUGAU-UUUACUUGGAUUCCCCUUCAGAACAUUUA
AF228365.1/1-73    CGGCAGAUGAUGAU-UUUACUUGGAUUCCCCUUCAGAACAUUU
#=GC SS_cons       <<<<_______..________.__._.______.___.___.___
//
```

The `#=GC SS_cons` line is the 2D structure consensus annotation, which represents the secondary structure of the RNA alignment. While not strictly required, it is recommended for Rfam alignments.

For more information, see the [Stockholm format specification](https://en.wikipedia.org/wiki/Stockholm_format).

## Contact us

If you have any questions or feedback, feel free to submit a GitHub issue or [email us](https://docs.rfam.org/en/latest/).
