# rfam-msa-qa
Tools for ensuring the quality of multiple sequence alignments (MSAs) in Rfam

## Stockholm File Validation

This repository includes a validation script for Stockholm format alignment files (`.so`, `.sto`, `.stk`).

### Usage

```bash
python3 validate_stockholm.py [-v] <file1.so> [file2.so ...]
```

Options:
- `-v, --verbose`: Print detailed validation information

### What is validated?

The script checks for:
- Required `# STOCKHOLM 1.0` header
- Required `//` terminator
- Presence of at least one sequence
- Consistent sequence lengths across the alignment

### Example

```bash
# Validate a single file
python3 validate_stockholm.py example_valid.so

# Validate multiple files with verbose output
python3 validate_stockholm.py -v file1.so file2.so file3.so
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
seq1    ACGU-ACGU-ACGU
seq2    ACGU-ACGU-ACGU
//
```

For more information, see the [Stockholm format specification](https://en.wikipedia.org/wiki/Stockholm_format).
