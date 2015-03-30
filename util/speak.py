"""
uitl.speak : A set of functions to handle the verbal output.

This module provides functions to help create or suppress verbal output. 

Functions list:
    * notquiet( args, words )
    * verbose( args, words )
"""


def notquiet( args, words ) :
    """Print statment if no quiet option set.

    This will print a statement as long as the quiet option was not passed.

    Parameters
    ----------
    args : list
        List of parsed arguments. args.quiet should be a bool. 
    words : str
        The string to print.

    Examples
    --------
    >>> args.quiet = True
    >>> notquiet(args, "Nothing should happen.")

    >>> args.quiet = False
    >>> notquiet(args, "Will print this.")
    Will print this.
    """
    if not args.quiet :
        print( words )

def verbose( args, words ) :
    """Print statment if verbose option set.

    This will print a statement when the verbose option was passed.

    Parameters
    ----------
    args : list
        List of parsed arguments. args.verbose should be a bool. 
    words : str
        The string to print.

    Examples
    --------
    >>> args.verbose = True
    >>> notquiet(args, "Will print this.")
    Will print this.

    >>> args.verbose = False
    >>> notquiet(args, "But not this.")
    
    """
    if args.verbose :
        print( words )
