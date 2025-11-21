#!/usr/bin/env python3
"""
Simple Stockholm format validator for MSA files.

This script validates Stockholm format alignment files (.so, .sto, .stk)
by checking for required structural elements and basic formatting.
"""

import sys
import argparse
from pathlib import Path


def validate_stockholm_file(filepath):
    """
    Validate a Stockholm format alignment file.
    
    Args:
        filepath: Path to the Stockholm file to validate
        
    Returns:
        tuple: (bool, list) - (is_valid, list_of_errors)
    """
    errors = []
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        return False, [f"Failed to read file: {e}"]
    
    if not lines:
        return False, ["File is empty"]
    
    # Check for Stockholm header
    header_found = False
    for i, line in enumerate(lines):
        if line.strip().startswith('# STOCKHOLM'):
            header_found = True
            # Validate version (typically 1.0)
            parts = line.strip().split()
            if len(parts) >= 3:
                version = parts[2]
                if version != '1.0':
                    errors.append(f"Line {i+1}: Non-standard Stockholm version '{version}' (expected 1.0)")
            break
    
    if not header_found:
        errors.append("Missing required '# STOCKHOLM 1.0' header")
    
    # Check for terminator
    terminator_found = False
    for i, line in enumerate(lines):
        if line.strip() == '//':
            terminator_found = True
            break
    
    if not terminator_found:
        errors.append("Missing required '//' terminator")
    
    # Check for at least one sequence
    sequence_count = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        # Skip empty lines, comments, markup, and terminator
        if not stripped or stripped.startswith('#') or stripped == '//':
            continue
        # Sequence lines should have at least two fields (name and sequence)
        parts = stripped.split()
        if len(parts) >= 2:
            sequence_count += 1
    
    if sequence_count == 0:
        errors.append("No sequences found in alignment")
    
    # Validate alignment consistency (all sequences should have same length)
    sequences = {}
    for i, line in enumerate(lines):
        stripped = line.strip()
        # Skip empty lines, comments, markup, and terminator
        if not stripped or stripped.startswith('#') or stripped == '//':
            continue
        
        parts = stripped.split(None, 1)  # Split on first whitespace only
        if len(parts) >= 2:
            seq_name = parts[0]
            seq_data = parts[1].replace(' ', '')  # Remove internal spaces
            
            if seq_name in sequences:
                # Append to existing sequence (Stockholm can have interleaved format)
                sequences[seq_name] += seq_data
            else:
                sequences[seq_name] = seq_data
    
    # Check all sequences have the same length
    if sequences:
        lengths = {name: len(seq) for name, seq in sequences.items()}
        unique_lengths = set(lengths.values())
        
        if len(unique_lengths) > 1:
            errors.append(f"Sequences have inconsistent lengths: {lengths}")
    
    is_valid = len(errors) == 0
    return is_valid, errors


def main():
    """Main function to run validation from command line."""
    parser = argparse.ArgumentParser(
        description='Validate Stockholm format MSA files'
    )
    parser.add_argument(
        'files',
        nargs='+',
        help='Stockholm file(s) to validate'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print detailed validation information'
    )
    
    args = parser.parse_args()
    
    all_valid = True
    
    for filepath in args.files:
        path = Path(filepath)
        
        if not path.exists():
            print(f"ERROR: File not found: {filepath}")
            all_valid = False
            continue
        
        is_valid, errors = validate_stockholm_file(filepath)
        
        if is_valid:
            if args.verbose:
                print(f"✓ {filepath}: Valid Stockholm file")
        else:
            print(f"✗ {filepath}: INVALID")
            for error in errors:
                print(f"  - {error}")
            all_valid = False
    
    if all_valid:
        print(f"\nAll {len(args.files)} file(s) passed validation.")
        return 0
    else:
        print(f"\nValidation failed for one or more files.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
