def compareHits(hitA, hitB):
    '''Returns:
    -   1 if hitB is a better hit than hitA for the purposes of
    taxonomic classification.
    -   -1 if hitA is better than hitB.
    -   0 if both hits identical
    '''
    if hitA.percent_identity < hitB.percent_identity:
        return 1
    elif hitA.percent_identity > hitB.percent_identity:
        return -1
    else:
        return 0