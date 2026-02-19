"""
Warning validators for Stockholm files.

These are non-critical issues that don't prevent validation.
"""


def check_ss_cons_from_parsed(gc_annotations, sequence_entries):
    """
    Check for 2D structure consensus annotation using pre-parsed merged data.

    Args:
        gc_annotations: dict of {feature: merged_data} from parser
        sequence_entries: list of (seq_name, full_seq_data) tuples (merged)

    Returns:
        tuple: (has_warning, warning_message_or_list)
    """
    ss_cons_data = gc_annotations.get('SS_cons')

    if ss_cons_data is None:
        return True, "Missing 2D structure consensus annotation (#=GC SS_cons)"

    warnings = []

    # Check length matches sequences
    ss_len = len(ss_cons_data)
    if sequence_entries:
        expected_len = len(sequence_entries[0][1])
        if ss_len != expected_len:
            warnings.append(
                f"#=GC SS_cons length ({ss_len}) does not match "
                f"sequence length ({expected_len})"
            )

    # Validate WUSS notation characters and bracket balancing
    ss_warnings = validate_ss_cons_format(ss_cons_data)
    warnings.extend(ss_warnings)

    if warnings:
        return True, warnings
    return False, None


def check_ss_cons(lines):
    """
    Check for 2D structure consensus annotation: existence and length.

    Args:
        lines: List of file lines

    Returns:
        tuple: (has_warning, warning_message_or_list)
            Returns a single string or a list of strings for multiple warnings.
    """
    ss_cons_data = None
    seq_lengths = []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith('#=GC SS_cons'):
            parts = stripped.split(None, 2)
            if len(parts) >= 3:
                ss_cons_data = parts[2]
        elif not stripped.startswith('#') and stripped and stripped != '//':
            parts = stripped.split(None, 1)
            if len(parts) >= 2:
                seq_lengths.append(len(parts[1].replace(' ', '')))

    if ss_cons_data is None:
        return True, "Missing 2D structure consensus annotation (#=GC SS_cons)"

    warnings = []

    # Check length matches sequences
    ss_len = len(ss_cons_data)
    if seq_lengths:
        expected_len = seq_lengths[0]
        if ss_len != expected_len:
            warnings.append(
                f"#=GC SS_cons length ({ss_len}) does not match "
                f"sequence length ({expected_len})"
            )

    # Validate WUSS notation characters and bracket balancing
    ss_warnings = validate_ss_cons_format(ss_cons_data)
    warnings.extend(ss_warnings)

    if warnings:
        return True, warnings
    return False, None


# Valid WUSS (Washington University Secondary Structure) characters
# Paired: () <> [] {} Aa Bb Cc Dd
# Unpaired: . , : _ - ~
WUSS_VALID_CHARS = set('.<>()[]{}AaBbCcDd,:_-~')

# Bracket pairs for balancing checks
WUSS_BRACKET_PAIRS = {
    '(': ')', '<': '>', '[': ']', '{': '}',
    'A': 'a', 'B': 'b', 'C': 'c', 'D': 'd',
}
WUSS_CLOSE_TO_OPEN = {v: k for k, v in WUSS_BRACKET_PAIRS.items()}


def validate_ss_cons_format(ss_cons_data):
    """
    Validate that #=GC SS_cons uses valid WUSS notation and has balanced brackets.

    Args:
        ss_cons_data: The SS_cons annotation string

    Returns:
        list: Warning messages (empty if valid)
    """
    warnings = []

    # Check for invalid characters
    invalid_chars = {}
    for i, ch in enumerate(ss_cons_data):
        if ch not in WUSS_VALID_CHARS:
            if ch not in invalid_chars:
                invalid_chars[ch] = i + 1  # 1-based position

    if invalid_chars:
        chars_str = ', '.join(
            f"'{ch}' (position {pos})" for ch, pos in invalid_chars.items()
        )
        warnings.append(f"#=GC SS_cons contains invalid characters: {chars_str}")

    # Check bracket balancing
    stacks = {}  # one stack per opening bracket type
    for opener in WUSS_BRACKET_PAIRS:
        stacks[opener] = []

    for i, ch in enumerate(ss_cons_data):
        if ch in WUSS_BRACKET_PAIRS:
            stacks[ch].append(i + 1)
        elif ch in WUSS_CLOSE_TO_OPEN:
            opener = WUSS_CLOSE_TO_OPEN[ch]
            if stacks[opener]:
                stacks[opener].pop()
            else:
                warnings.append(
                    f"#=GC SS_cons has unmatched '{ch}' at position {i + 1}"
                )

    for opener, stack in stacks.items():
        if stack:
            closer = WUSS_BRACKET_PAIRS[opener]
            warnings.append(
                f"#=GC SS_cons has {len(stack)} unmatched '{opener}' "
                f"(missing '{closer}')"
            )

    return warnings


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
