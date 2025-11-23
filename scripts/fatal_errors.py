"""
Fatal error validators for Stockholm files.

These errors cannot be automatically fixed and must terminate execution.
"""


def check_header(lines):
    """
    Check for required Stockholm header.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (is_valid, error_message, line_number)
    """
    for i, line in enumerate(lines):
        if line.strip().startswith('# STOCKHOLM'):
            parts = line.strip().split()
            if len(parts) >= 3:
                version = parts[2]
                if version != '1.0':
                    return False, f"Non-standard Stockholm version '{version}' (expected 1.0)", i + 1
            return True, None, None
    return False, "Missing required '# STOCKHOLM 1.0' header", None


def check_terminator(lines):
    """
    Check for required terminator.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (is_valid, error_message)
    """
    for line in lines:
        if line.strip() == '//':
            return True, None
    return False, "Missing required '//' terminator"


def check_empty_file(lines):
    """
    Check if file is empty.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not lines:
        return False, "File is empty"
    return True, None


def check_no_sequences(sequences):
    """
    Check if alignment has sequences.
    
    Args:
        sequences: Dictionary of sequences
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not sequences:
        return False, "No sequences found in alignment"
    return True, None


def check_sequence_lengths(sequences):
    """
    Check that all sequences have the same length.
    
    Args:
        sequences: Dictionary of {seq_name: seq_data}
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not sequences:
        return True, None
    
    lengths = {name: len(seq) for name, seq in sequences.items()}
    unique_lengths = set(lengths.values())
    
    if len(unique_lengths) > 1:
        return False, f"Sequences have inconsistent lengths: {lengths}"
    return True, None


def check_sequence_characters(sequences):
    """
    Check that sequences contain only valid characters (no whitespace).
    
    Args:
        sequences: Dictionary of {seq_name: seq_data}
        
    Returns:
        tuple: (is_valid, error_message)
    """
    for seq_name, seq_data in sequences.items():
        # Check for whitespace
        if any(c.isspace() for c in seq_data):
            return False, f"Sequence '{seq_name}' contains whitespace characters"
    return True, None


def validate_gap_characters(sequences):
    """
    Validate that gaps are indicated by . or -.
    
    This is informational - we accept gaps as either . or -
    
    Args:
        sequences: Dictionary of {seq_name: seq_data}
        
    Returns:
        tuple: (is_valid, info_message)
    """
    # This is not a fatal error, just validation
    # Gaps can be . or - and both are valid
    return True, None
