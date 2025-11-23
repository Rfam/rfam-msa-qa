"""
Fixable error handlers for Stockholm files.

These errors can be automatically corrected.
"""


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


def find_duplicates_from_entries(sequence_entries):
    """
    Find duplicate sequences from a list of sequence entries.
    
    Args:
        sequence_entries: List of (seq_name, seq_data) tuples
        
    Returns:
        tuple: (unique_sequences_dict, duplicate_indices)
    """
    seen_keys = {}
    unique_sequences = {}
    duplicate_indices = []
    
    for idx, (seq_name, seq_data) in enumerate(sequence_entries):
        accession, coords = parse_sequence_identifier(seq_name)
        key = (accession, coords, seq_data)
        
        if key not in seen_keys:
            seen_keys[key] = idx
            unique_sequences[seq_name] = seq_data
        else:
            duplicate_indices.append(idx)
    
    return unique_sequences, duplicate_indices


def remove_duplicates(lines):
    """
    Remove duplicate sequences from Stockholm file lines.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (corrected_lines, num_duplicates_removed)
    """
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
    
    # Find duplicates
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
    
    # Reconstruct file
    corrected_lines = []
    
    # Add header lines
    corrected_lines.extend(header_lines)
    
    # Calculate max sequence name length for formatting
    max_name_len = max(len(name) for name, _ in unique_entries) if unique_entries else 30
    name_width = max(max_name_len + 2, 30)
    
    # Add unique sequences
    for seq_name, seq_content in unique_entries:
        corrected_lines.append(f"{seq_name:<{name_width}} {seq_content}\n")
    
    # Add annotation lines
    corrected_lines.extend(annotation_lines)
    
    # Add footer lines
    corrected_lines.extend(footer_lines)
    
    return corrected_lines, duplicates_count
