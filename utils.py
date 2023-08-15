
from math import gcd
import typing


def reduce_seq_to_lowest(seq: typing.Sequence[int]):
    """Returns the smallest ratio of the given sequence of integers"""
    g = gcd(*seq)
    return (i / g for i in seq)

