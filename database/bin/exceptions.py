#!/usr/bin/env python3
"""
Custom exceptions for CG2AT.

Replaces bare sys.exit() calls so errors can be caught, tested, and
re-raised cleanly by callers rather than killing the interpreter mid-call.

Usage
-----
    from exceptions import CG2ATError, FragmentNotFoundError, \
        TopologyError, InputError

    raise FragmentNotFoundError(f"Cannot find fragment: {residue}/{name}.pdb")

In main() / entry-points wrap with:

    try:
        run()
    except CG2ATError as exc:
        sys.exit(str(exc))
"""


class CG2ATError(RuntimeError):
    """Base exception for all CG2AT errors."""


class InputError(CG2ATError):
    """Raised when user-supplied input files or flags are invalid."""


class FragmentNotFoundError(CG2ATError):
    """Raised when a residue fragment cannot be located in the database."""


class TopologyError(CG2ATError):
    """Raised when a topology or force-field file is malformed or missing."""


class GromacsError(CG2ATError):
    """Raised when a GROMACS sub-process fails or cannot be found."""


class SequenceError(CG2ATError):
    """Raised when CG and atomistic sequences cannot be reconciled."""
