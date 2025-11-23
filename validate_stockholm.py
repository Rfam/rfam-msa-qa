#!/usr/bin/env python3
"""
Simple Stockholm format validator for MSA files.

This script validates Stockholm format alignment files (.so, .sto, .stk)
by checking for required structural elements and basic formatting.
It also removes duplicate sequences from the alignment.
"""

import sys
import argparse
from pathlib import Path


def parse_sequence_identifier(seq_name):
    """
    Parse sequence identifier to extract accession and coordinates.
    
    Format: ACCESSION/START-END (e.g., AF228364.1/1-74)
    
    Args:
        seq_name: Sequence identifier string
        
    Returns:
        tuple: (accession, coordinates) or (seq_name, None) if no coordinates found
    """
    if '/' in seq_name:
        parts = seq_name.split('/', 1)
        return parts[0], parts[1]
    return seq_name, None


def validate_stockholm_file(filepath, check_duplicates=False):
    """
    Validate a Stockholm format alignment file.
    
    Args:
        filepath: Path to the Stockholm file to validate
        check_duplicates: If True, count duplicate sequences
        
    Returns:
        tuple: (bool, list, int, list) - (is_valid, list_of_errors, num_duplicates_found, list_of_warnings)
    """
    errors = []
    warnings = []
    num_duplicates = 0
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        return False, [f"Failed to read file: {e}"], 0, []
    
    if not lines:
        return False, ["File is empty"], 0, []
    
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
    
    # Check for 2D structure consensus annotation
    ss_cons_found = False
    for line in lines:
        if line.strip().startswith('#=GC SS_cons'):
            ss_cons_found = True
            break
    
    if not ss_cons_found:
        warnings.append("Missing 2D structure consensus annotation (#=GC SS_cons)")
    
    # Parse sequences
    sequences = {}
    for i, line in enumerate(lines):
        stripped = line.strip()
        # Skip empty lines, comments, markup, and terminator
        if not stripped or stripped.startswith('#') or stripped == '//':
            continue
        
        parts = stripped.split(None, 1)  # Split on first whitespace only
        if len(parts) >= 2:
            seq_name = parts[0]
            # Remove internal spaces - Stockholm format allows spaces in sequence data
            seq_data = parts[1].replace(' ', '')
            
            if seq_name in sequences:
                # Append to existing sequence (Stockholm can have interleaved format)
                sequences[seq_name] += seq_data
            else:
                sequences[seq_name] = seq_data
    
    if not sequences:
        errors.append("No sequences found in alignment")
    
    # Check for duplicates if requested
    if check_duplicates and sequences:
        unique_sequences = {}
        duplicates_found = []
        
        for seq_name, seq_data in sequences.items():
            accession, coords = parse_sequence_identifier(seq_name)
            key = (accession, coords, seq_data)
            
            if key not in unique_sequences:
                unique_sequences[key] = seq_name
            else:
                duplicates_found.append(seq_name)
        
        num_duplicates = len(duplicates_found)
    
    is_valid = len(errors) == 0
    return is_valid, errors, num_duplicates, warnings


def remove_duplicates_from_file(filepath, output_filepath=None):
    """
    Remove duplicate sequences from a Stockholm file.
    
    Duplicates are defined as sequences with the same accession/identifier,
    coordinates, and exact sequence data.
    
    Args:
        filepath: Path to input Stockholm file
        output_filepath: Path to output file (if None, overwrites input)
        
    Returns:
        int: Number of duplicates removed
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading file: {e}")
        return 0
    
    # Separate different types of lines
    header_lines = []
    annotation_lines = []
    footer_lines = []
    sequence_entries = []  # List of (seq_name, seq_content) tuples
    
    sequences_started = False
    sequences_ended = False
    
    for line in lines:
        stripped = line.strip()
        
        if stripped == '//':
            sequences_ended = True
            footer_lines.append(line)
        elif sequences_ended:
            footer_lines.append(line)
        elif stripped.startswith('# STOCKHOLM'):
            header_lines.append(line)
            sequences_started = True
        elif stripped.startswith('#=GF') or stripped.startswith('#=GC') or stripped.startswith('#=GS') or stripped.startswith('#=GR'):
            annotation_lines.append(line)
        elif stripped.startswith('#'):
            if not sequences_started:
                header_lines.append(line)
            else:
                annotation_lines.append(line)
        elif stripped:
            # This is a sequence line
            parts = stripped.split(None, 1)
            if len(parts) >= 2:
                seq_name = parts[0]
                seq_content = parts[1].replace(' ', '')
                sequence_entries.append((seq_name, seq_content))
        else:
            # Empty line
            if not sequences_started:
                header_lines.append(line)
            elif not sequences_ended:
                # Keep empty lines in body
                pass
            else:
                footer_lines.append(line)
    
    # Find duplicates based on accession, coordinates, and sequence
    seen_keys = {}
    unique_entries = []
    duplicates_count = 0
    
    for seq_name, seq_content in sequence_entries:
        accession, coords = parse_sequence_identifier(seq_name)
        key = (accession, coords, seq_content)
        
        if key not in seen_keys:
            seen_keys[key] = True
            unique_entries.append((seq_name, seq_content))
        else:
            duplicates_count += 1
    
    # Write output file
    if output_filepath is None:
        output_filepath = filepath
    
    # Calculate max sequence name length for formatting
    max_name_len = max(len(name) for name, _ in unique_entries) if unique_entries else 30
    # Add some padding
    name_width = max(max_name_len + 2, 30)
    
    with open(output_filepath, 'w') as f:
        # Write header lines
        for line in header_lines:
            f.write(line)
        
        # Write unique sequences
        for seq_name, seq_content in unique_entries:
            f.write(f"{seq_name:<{name_width}} {seq_content}\n")
        
        # Write annotation lines
        for line in annotation_lines:
            f.write(line)
        
        # Write footer lines
        for line in footer_lines:
            f.write(line)
    
    return duplicates_count


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
    parser.add_argument(
        '--remove-duplicates',
        action='store_true',
        help='Remove duplicate sequences (same accession, coordinates, and sequence)'
    )
    
    args = parser.parse_args()
    
    all_valid = True
    
    for filepath in args.files:
        path = Path(filepath)
        
        if not path.exists():
            print(f"ERROR: File not found: {filepath}")
            all_valid = False
            continue
        
        # Remove duplicates if requested
        if args.remove_duplicates:
            num_removed = remove_duplicates_from_file(filepath)
            if num_removed > 0:
                print(f"Removed {num_removed} duplicate sequence(s) from {filepath}")
        
        is_valid, errors, num_duplicates, warnings = validate_stockholm_file(filepath, check_duplicates=not args.remove_duplicates)
        
        if is_valid:
            if args.verbose:
                print(f"✓ {filepath}: Valid Stockholm file")
                if not args.remove_duplicates and num_duplicates > 0:
                    print(f"  Note: Found {num_duplicates} duplicate sequence(s)")
            # Display warnings even if file is valid
            if warnings:
                if not args.verbose:
                    print(f"✓ {filepath}: Valid Stockholm file")
                for warning in warnings:
                    print(f"  ⚠ Warning: {warning}")
        else:
            print(f"✗ {filepath}: INVALID")
            for error in errors:
                print(f"  - {error}")
            # Display warnings even if there are errors
            if warnings:
                for warning in warnings:
                    print(f"  ⚠ Warning: {warning}")
            all_valid = False
    
    if all_valid:
        print(f"\nAll {len(args.files)} file(s) passed validation.")
        return 0
    else:
        print(f"\nValidation failed for one or more files.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
