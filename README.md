# rfam-msa-qa
Tools for ensuring the quality of multiple sequence alignments (MSAs) in Rfam

## Stockholm File Validation

This repository includes a validation script for Stockholm format alignment files (`.so`, `.sto`, `.stk`).

### Usage

```bash
python3 validate_stockholm.py [-v] [--remove-duplicates] <file1.so> [file2.so ...]
```

Options:
- `-v, --verbose`: Print detailed validation information
- `--remove-duplicates`: Remove duplicate sequences from the file(s)

### What is validated?

The script checks for:
- Required `# STOCKHOLM 1.0` header
- Required `//` terminator
- Presence of at least one sequence

The script also provides warnings for:
- Missing 2D structure consensus annotation (`#=GC SS_cons`)

Note: Sequences are allowed to have different lengths in the alignment.

### Sequence Format

Sequences should follow the format: `ACCESSION/START-END` where:
- `ACCESSION` is the sequence identifier (e.g., from GenBank like `AF228364.1`)
- `START-END` are the coordinates indicating which portion of the original sequence is included (e.g., `1-74`)

Example: `AF228364.1/1-74`

### Duplicate Detection and Removal

The script can detect and remove duplicate sequences using the `--remove-duplicates` flag. Duplicates are defined as sequences that have:
1. The same accession/identifier
2. The same coordinates
3. The exact same sequence data

### Example

```bash
# Validate a single file
python3 validate_stockholm.py example_valid.so

# Validate multiple files with verbose output
python3 validate_stockholm.py -v file1.so file2.so file3.so

# Remove duplicates from a file
python3 validate_stockholm.py --remove-duplicates file.so
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
