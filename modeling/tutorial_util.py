import numpy


def classify_crosslinks_by_score(xlres, number_of_classes):
    '''Creates the requested number of classes and partitions crosslinks
       according to their identification scores. Classes are defined in
       the psi key.'''
    if xlres.id_score_key is not None:
        scores=xlres.get_values(xlres.id_score_key)
    else:
        raise ValueError('The crosslink database does not contain score values')
    minscore=min(scores)
    maxscore=max(scores)
    scoreclasses=numpy.linspace(minscore, maxscore, number_of_classes+1)
    if xlres.psi_key is None:
        xlres.create_new_keyword(xlres.psi_key,values_from_keyword=None)
    for xl in xlres:
        score=xl[xlres.id_score_key]
        for n,classmin in enumerate(scoreclasses[0:-1]):
            if score>=classmin and score<=scoreclasses[n+1]:
                xl[xlres.psi_key]=str("CLASS_"+str(n))
    xlres._update()
