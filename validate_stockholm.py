#!/usr/bin/env python3
"""
Stockholm format validator for MSA files.

This script validates Stockholm format alignment files (.so, .sto, .stk)
by checking for required structural elements and basic formatting.
It can automatically fix certain errors like duplicate sequences.
"""

import sys
import argparse
from pathlib import Path

# Import validation modules
from scripts import fatal_errors, fixable_errors, warnings, parser


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


def fix_file(filepath, output_mode='file'):
    """
    Fix fixable errors in a Stockholm file.
    
    Args:
        filepath: Path to input Stockholm file
        output_mode: 'stdout' to print to stdout, 'file' to create _corrected file
        
    Returns:
        tuple: (success, num_fixes_applied, output_path)
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return False, 0, None
    
    # Apply fixes (currently only duplicate removal)
    corrected_lines, num_duplicates = fixable_errors.remove_duplicates(lines)
    
    if output_mode == 'stdout':
        # Print to stdout
        for line in corrected_lines:
            print(line, end='')
        return True, num_duplicates, None
    else:
        # Create _corrected file
        path = Path(filepath)
        output_path = path.parent / f"{path.stem}_corrected{path.suffix}"
        
        with open(output_path, 'w') as f:
            f.writelines(corrected_lines)
        
        return True, num_duplicates, str(output_path)


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
        
        # Handle fixing if requested
        if args.fix and result['can_be_fixed']:
            success, num_fixes, output_path = fix_file(filepath, args.output_mode)
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
