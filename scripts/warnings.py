"""
Warning validators for Stockholm files.

These are non-critical issues that don't prevent validation.
"""


def check_ss_cons(lines):
    """
    Check for 2D structure consensus annotation.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (has_warning, warning_message)
    """
    for line in lines:
        if line.strip().startswith('#=GC SS_cons'):
            return False, None
    return True, "Missing 2D structure consensus annotation (#=GC SS_cons)"


def check_line_length(lines, max_length=10000):
    """
    Check for lines exceeding Rfam alignment limits.
    
    Args:
        lines: List of file lines
        max_length: Maximum allowed line length
        
    Returns:
        tuple: (has_warning, warning_messages_list)
    """
    warnings = []
    for i, line in enumerate(lines):
        if len(line) > max_length:
            warnings.append(f"Line {i+1} exceeds {max_length} character limit (length: {len(line)})")
    
    if warnings:
        return True, warnings
    return False, []
