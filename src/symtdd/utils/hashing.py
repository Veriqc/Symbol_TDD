import math


def hash_double(d: float) -> int:
    """
    https://book.huihoo.com/data-structures-and-algorithms-with-object-oriented-design-patterns-in-c++/html/page217.html
    """
    return (
        math.floor((2 * math.fabs(math.frexp(d)[0]) - 1) * (2**64 - 1))
        if d != 0
        else 0
    )
