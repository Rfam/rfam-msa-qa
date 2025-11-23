"""
Utility functions for parsing Stockholm files.
"""


def parse_stockholm_file(lines):
    """
    Parse Stockholm file into structured components.
    
    This parser treats each sequence line individually and does NOT combine
    sequences with the same name (interleaved format). Instead, it keeps
    each occurrence separate for duplicate detection.
    
    Args:
        lines: List of file lines
        
    Returns:
        dict: Parsed file structure with sequences, headers, annotations, etc.
    """
    header_found = False
    terminator_found = False
    ss_cons_found = False
    sequence_entries = []  # List of (seq_name, seq_data) tuples
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        # Check for Stockholm header
        if stripped.startswith('# STOCKHOLM'):
            header_found = True
        
        # Check for terminator
        elif stripped == '//':
            terminator_found = True
        
        # Check for 2D structure consensus annotation
        elif stripped.startswith('#=GC SS_cons'):
            ss_cons_found = True
        
        # Parse sequence lines (skip empty, comments, markup, and terminator)
        elif stripped and not stripped.startswith('#'):
            parts = stripped.split(None, 1)  # Split on first whitespace only
            if len(parts) >= 2:
                seq_name = parts[0]
                # Remove internal spaces - Stockholm format allows spaces in sequence data
                seq_data = parts[1].replace(' ', '')
                sequence_entries.append((seq_name, seq_data))
    
    return {
        'header_found': header_found,
        'terminator_found': terminator_found,
        'ss_cons_found': ss_cons_found,
        'sequence_entries': sequence_entries,
        'lines': lines
    }
