from typing import Any, Dict


def compareHits(hitA: Dict[str, Any], hitB: Dict[str, Any]):
    '''Returns:
    -   1 if hitB is a better hit than hitA for the purposes of
    taxonomic classification.
    -   -1 if hitA is better than hitB.
    -   0 if both hits identical
    '''
    a = hitA.get('percent_identity', None)
    b = hitB.get('percent_identity', None)
    assert not a is None and not b is None
    if a < b:
        return 1
    elif a > b:
        return -1
    else:
        return 0