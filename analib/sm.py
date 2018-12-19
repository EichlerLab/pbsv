"""
Tools for managing Snakemake resources.
"""


def named_list_set(named_list, key, value, wildcards=None):
    """
    Set a value on a named list in Snakemake. This includes objects for wildcards, input, output, parameters, and log.

    :param named_list: Named list (wildcards, input, output, parameters, etc).
    :param key: Key.
    :param value: Value to set. May be a string or an input function taking a single argument, `wildcards`.
        This parameter must be set if `value` is a function.
    :param wildcards: Format `value` with `wildcards` if set. If `value` is a function, call it with
        `wildcards` as its only parameter.
    """

    if callable(value):
        # Input function
        if wildcards is None:
            raise RuntimeError('Cannot execute input function with wildcards = None')

        value = value(wildcards)

    else:
        # Format
        if wildcards is not None:
            value = value.format(**wildcards)

    # Add key if missing
    if key not in named_list.keys().keys():
        named_list.append(None)
        named_list.add_name(key)

    # Set value
    setattr(named_list, key, value)
    named_list[named_list._names[key][0]] = value
