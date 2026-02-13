#!/usr/bin/env python3
"""
Stockholm format validator for MSA files.

This script validates Stockholm format alignment files (.so, .sto, .stk)
by checking for required structural elements and basic formatting.
It can automatically fix certain errors like duplicate sequences.
"""

# Suppress BioPython warnings before any imports
import warnings as stdlib_warnings
stdlib_warnings.filterwarnings("ignore", module="Bio")

import sys
import argparse
from pathlib import Path

# Import validation modules
from scripts import fatal_errors, fixable_errors, stockholm_warnings as warnings, parser


def validate_stockholm_file(filepath):
    """
    Validate a Stockholm format alignment file.
    
    Args:
        filepath: Path to the Stockholm file to validate
        
    Returns:
        dict: Validation results with errors, warnings, and fixable issues
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        return {
            'is_valid': False,
            'fatal_errors': [f"Failed to read file: {e}"],
            'fixable_errors': [],
            'warnings': [],
            'can_be_fixed': False
        }
    
    # Parse the file
    parsed = parser.parse_stockholm_file(lines)
    sequence_entries = parsed['sequence_entries']
    
    # Build unique sequences dict for validation
    unique_sequences, duplicates = fixable_errors.find_duplicates_from_entries(sequence_entries)
    
    # Collect errors and warnings
    fatal_error_list = []
    fixable_error_list = []
    warning_list = []
    
    # Check for fatal errors
    is_valid, error_msg = fatal_errors.check_empty_file(lines)
    if not is_valid:
        fatal_error_list.append(error_msg)
        return {
            'is_valid': False,
            'fatal_errors': fatal_error_list,
            'fixable_errors': [],
            'warnings': [],
            'can_be_fixed': False
        }
    
    is_valid, error_msg, line_num = fatal_errors.check_header(lines)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    is_valid, error_msg = fatal_errors.check_terminator(lines)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    is_valid, error_msg = fatal_errors.check_no_sequences(unique_sequences)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    # Check sequence length consistency
    is_valid, error_msg = fatal_errors.check_sequence_lengths(unique_sequences)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    # Check sequence characters
    is_valid, error_msg = fatal_errors.check_sequence_characters(unique_sequences)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    # Check for fixable errors
    if len(duplicates) > 0:
        fixable_error_list.append(f"Found {len(duplicates)} duplicate sequence(s)")

    # Check for missing coordinates
    missing_coords = fixable_errors.find_missing_coordinates(sequence_entries)
    if missing_coords:
        fixable_error_list.append(f"Found {len(missing_coords)} sequence(s) missing coordinates")

    # Check for overlapping sequences
    overlapping = fixable_errors.find_overlapping_sequences(sequence_entries)
    if overlapping:
        fixable_error_list.append(f"Found {len(overlapping)} overlapping sequence(s)")

    # Check for warnings
    has_warning, warning_msg = warnings.check_ss_cons(lines)
    if has_warning:
        warning_list.append(warning_msg)
    
    has_warning, warning_msgs = warnings.check_line_length(lines)
    if has_warning:
        warning_list.extend(warning_msgs)
    
    # Determine if file can be fixed
    can_be_fixed = len(fatal_error_list) == 0 and len(fixable_error_list) > 0
    
    return {
        'is_valid': len(fatal_error_list) == 0,
        'fatal_errors': fatal_error_list,
        'fixable_errors': fixable_error_list,
        'warnings': warning_list,
        'can_be_fixed': can_be_fixed,
        'lines': lines
    }


def fix_file(filepath, output_mode='file', verbose=False):
    """
    Fix fixable errors in a Stockholm file.

    Args:
        filepath: Path to input Stockholm file
        output_mode: 'stdout' to print to stdout, 'file' to create _corrected file
        verbose: Print progress messages

    Returns:
        tuple: (success, num_fixes_applied, output_path)
    """
    import re

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return False, 0, None

    # Parse to check for missing coordinates
    parsed = parser.parse_stockholm_file(lines)
    sequence_entries = parsed['sequence_entries']
    missing_coords = fixable_errors.find_missing_coordinates(sequence_entries)

    num_fixes = 0
    id_mapping = {}
    to_remove = set()

    # Fix missing coordinates if any
    if missing_coords:
        if verbose:
            print(f"  Fixing {len(missing_coords)} missing coordinates (downloading from NCBI)...")
        id_mapping, to_remove = fixable_errors.fix_missing_coordinates(filepath, verbose=verbose)
        # Count sequences that got new coordinates (where key != value)
        num_coords_fixed = len([k for k, v in id_mapping.items() if k != v])
        num_fixes += num_coords_fixed
        if verbose and to_remove:
            print(f"  {len(to_remove)} sequence(s) will be removed (no accurate match)")

    # Apply coordinate fixes, remove bad sequences, and remove duplicates
    header_lines = []
    gf_lines = []  # #=GF lines go at top after header
    gs_gr_lines = []  # #=GS/#=GR lines (sequence annotations)
    gc_lines = []  # #=GC lines go at bottom before //
    footer_lines = []
    processed_sequences = []  # List of (new_name, seq_content)
    seen_sequences = {}
    num_duplicates = 0
    num_removed = 0
    in_footer = False

    for line in lines:
        stripped = line.strip()

        if stripped == '//':
            in_footer = True
            footer_lines.append(line)
        elif in_footer:
            footer_lines.append(line)
        elif stripped.startswith('# STOCKHOLM'):
            header_lines.append(line)
        elif stripped.startswith('#=GF'):
            # File-level annotations - keep at top
            gf_lines.append(line)
        elif stripped.startswith('#=GC'):
            # Column-level annotations - keep at bottom
            gc_lines.append(line)
        elif stripped.startswith('#=GS') or stripped.startswith('#=GR'):
            # Sequence-level annotations - update accession if needed
            parts = re.split(r'\s+', stripped)
            if len(parts) >= 3:
                seq_id = parts[1]
                # Skip annotation for removed sequences
                if seq_id in to_remove:
                    continue
                # Update accession if mapped
                if seq_id in id_mapping:
                    parts[1] = id_mapping[seq_id]
                    info = ' '.join(parts[:-1])
                    gs_gr_lines.append((seq_id, f"{info:<30} {parts[-1]}\n"))
                else:
                    gs_gr_lines.append((seq_id, line))
            else:
                gs_gr_lines.append((None, line))
        elif stripped.startswith('#') or line == "\n":
            # Other comments or blank lines - keep at top
            header_lines.append(line)
        elif stripped:
            # Sequence line
            parts = stripped.split(None, 1)
            if len(parts) >= 2:
                seq_name = parts[0]
                seq_content = parts[1].replace(' ', '')

                # Skip sequences marked for removal
                if seq_name in to_remove:
                    num_removed += 1
                    continue

                # Apply coordinate mapping if available
                if seq_name in id_mapping:
                    new_name = id_mapping[seq_name]
                else:
                    new_name = seq_name

                # Check for duplicates
                accession, coords = fixable_errors.parse_sequence_identifier(new_name)
                key = (accession, coords, seq_content)

                if key not in seen_sequences:
                    seen_sequences[key] = True
                    processed_sequences.append((new_name, seq_content))
                else:
                    num_duplicates += 1

    num_fixes += num_duplicates

    # Check for overlapping sequences and remove them
    overlapping = fixable_errors.find_overlapping_sequences(processed_sequences)
    if overlapping:
        if verbose:
            print(f"  Found {len(overlapping)} overlapping sequence(s) to remove")
        # Filter out overlapping sequences
        processed_sequences = [(name, seq) for name, seq in processed_sequences if name not in overlapping]
        num_removed += len(overlapping)

    # Validate sequences against NCBI (with BLAST fallback)
    if verbose:
        print(f"  Validating {len(processed_sequences)} sequence(s) against NCBI...")
    invalid_seqs, not_found_seqs, mismatched_seqs, blast_fixed = fixable_errors.validate_sequences_against_ncbi(processed_sequences, verbose=verbose)

    # Apply BLAST fixes: update sequence names with BLAST-found accession/coords
    if blast_fixed:
        if verbose:
            print(f"  {len(blast_fixed)} sequence(s) rescued via BLAST")
        processed_sequences = [
            (blast_fixed[name], seq) if name in blast_fixed else (name, seq)
            for name, seq in processed_sequences
        ]
        num_fixes += len(blast_fixed)

    if invalid_seqs:
        if not_found_seqs:
            print(f"  {len(not_found_seqs)} sequence(s) not found in NCBI (invalid accession?)")
        if mismatched_seqs:
            print(f"  {len(mismatched_seqs)} sequence(s) don't match NCBI data")
        processed_sequences = [(name, seq) for name, seq in processed_sequences if name not in invalid_seqs]
        num_removed += len(invalid_seqs)

    # Build final corrected lines in proper Stockholm order
    corrected_lines = []

    # 1. Header (# STOCKHOLM 1.0)
    corrected_lines.extend(header_lines)

    # 2. File-level annotations (#=GF) - keep at top
    corrected_lines.extend(gf_lines)

    # 3. Sequence-level annotations (#=GS) for kept sequences
    kept_seq_names = {name for name, _ in processed_sequences}
    for orig_seq_id, ann_line in gs_gr_lines:
        if orig_seq_id is None:
            corrected_lines.append(ann_line)
        elif orig_seq_id in id_mapping:
            mapped_name = id_mapping[orig_seq_id]
            if mapped_name in kept_seq_names:
                corrected_lines.append(ann_line)
        elif orig_seq_id in kept_seq_names:
            corrected_lines.append(ann_line)

    # 4. Sequences
    for new_name, seq_content in processed_sequences:
        corrected_lines.append(f"{new_name:<30} {seq_content}\n")

    # 5. Column-level annotations (#=GC) - after sequences
    corrected_lines.extend(gc_lines)

    corrected_lines.extend(footer_lines)

    if verbose and num_removed > 0:
        print(f"  Removed {num_removed} sequence(s) total from output")

    if output_mode == 'stdout':
        for line in corrected_lines:
            print(line, end='')
        return True, num_fixes, None
    else:
        path = Path(filepath)
        output_path = path.parent / f"{path.stem}_corrected{path.suffix}"

        with open(output_path, 'w') as f:
            f.writelines(corrected_lines)

        return True, num_fixes, str(output_path)


def main():
    """Main function to run validation from command line."""
    parser_arg = argparse.ArgumentParser(
        description='Validate Stockholm format MSA files'
    )
    parser_arg.add_argument(
        'files',
        nargs='+',
        help='Stockholm file(s) to validate'
    )
    parser_arg.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print detailed validation information'
    )
    parser_arg.add_argument(
        '--fix',
        action='store_true',
        help='Attempt to fix fixable errors'
    )
    parser_arg.add_argument(
        '--output-mode',
        choices=['stdout', 'file'],
        default='file',
        help='Output mode for fixed files: stdout or create _corrected file (default: file)'
    )
    
    args = parser_arg.parse_args()
    
    all_valid = True
    
    for filepath in args.files:
        path = Path(filepath)
        
        if not path.exists():
            print(f"ERROR: File not found: {filepath}")
            all_valid = False
            continue
        
        # Validate the file
        result = validate_stockholm_file(filepath)
        
        # Handle fixing if requested (always run to validate sequences against NCBI)
        if args.fix and result['is_valid']:
            success, num_fixes, output_path = fix_file(filepath, args.output_mode, verbose=args.verbose)
            if success:
                if args.output_mode == 'file':
                    print(f"✓ Fixed {num_fixes} issue(s) in {filepath}")
                    print(f"  Corrected file saved to: {output_path}")
                else:
                    # Already printed to stdout
                    pass
                
                # Re-validate the fixed content
                if args.output_mode == 'file':
                    result = validate_stockholm_file(output_path)
        
        # Display results
        if result['is_valid'] and not result['fixable_errors']:
            if args.verbose or result['warnings']:
                print(f"✓ {filepath}: Valid Stockholm file")
            
            # Display warnings
            for warning in result['warnings']:
                print(f"  ⚠ Warning: {warning}")
        
        elif result['is_valid'] and result['fixable_errors']:
            print(f"⚠ {filepath}: Valid but has fixable issues")
            for error in result['fixable_errors']:
                print(f"  • {error}")
            
            # Display warnings
            for warning in result['warnings']:
                print(f"  ⚠ Warning: {warning}")
            
            if not args.fix:
                print(f"  Tip: Use --fix to automatically correct these issues")
        
        else:
            print(f"✗ {filepath}: INVALID")
            for error in result['fatal_errors']:
                print(f"  - {error}")
            
            # Display fixable errors
            for error in result['fixable_errors']:
                print(f"  • Fixable: {error}")
            
            # Display warnings
            for warning in result['warnings']:
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
