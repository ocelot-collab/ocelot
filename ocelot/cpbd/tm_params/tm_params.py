class TMParams:
    """
    Lightweight data container passed from an atom hook to a transformation.

    ``TMParams`` subclasses carry the precomputed element-specific data that a
    transformation needs at runtime. Some subclasses are pure data holders,
    while others add small helpers such as rotated-matrix accessors.
    """

    def __init__(self) -> None:
        pass
