import IMP.pmi.io.crosslink
import IMP.pmi.macros
import IMP.atom
import numpy


# Monkey patch CrossLinkDataBase.append_database to incorporate a post-2.11
# bug fix
def _append_database(self,CrossLinkDataBase2):
    name1=self.get_name()
    name2=CrossLinkDataBase2.get_name()
    if name1 == name2:
        name1=id(self)
        name2=id(CrossLinkDataBase2)
        self.set_name(name1)
        CrossLinkDataBase2.set_name(name2)
    #rename first database:
    new_data_base={}
    for k in self.data_base:
        new_data_base[k]=self.data_base[k]
    for k in CrossLinkDataBase2.data_base:
        new_data_base[k]=CrossLinkDataBase2.data_base[k]
    self.data_base=new_data_base
    self._update()


# Monkey patch AnalysisReplicaExchange.align to incorporate a post-2.11
# bug fix
def _align(self):
    tr = IMP.atom.get_transformation_aligning_first_to_second(self.sel1_alignment, self.sel0_alignment)

    for rb in self.rbs1:
        IMP.core.transform(rb, tr)
        for bead in self.beads1:
            try:
                IMP.core.transform(IMP.core.XYZ(bead), tr)
            except:
                continue
        self.model.update()


# This method isn't in IMP 2.11, so patch it in
def _classify_crosslinks_by_score(self, number_of_classes):
    '''Creates the requested number of classes and partitions crosslinks
       according to their identification scores. Classes are defined in
       the psi key.'''
    if self.id_score_key is not None:
        scores=self.get_values(self.id_score_key)
    else:
        raise ValueError('The crosslink database does not contain score values')
    minscore=min(scores)
    maxscore=max(scores)
    scoreclasses=numpy.linspace(minscore, maxscore, number_of_classes+1)
    if self.psi_key is None:
        self.create_new_keyword(self.psi_key,values_from_keyword=None)
    for xl in self:
        score=xl[self.id_score_key]
        for n,classmin in enumerate(scoreclasses[0:-1]):
            if score>=classmin and score<=scoreclasses[n+1]:
                xl[self.psi_key]=str("CLASS_"+str(n))
    self._update()


if IMP.__version__ == '2.11.1':
    IMP.pmi.macros.AnalysisReplicaExchange.align = _align
    IMP.pmi.io.crosslink.CrossLinkDataBase.append_database = _append_database
    IMP.pmi.io.crosslink.CrossLinkDataBase.classify_crosslinks_by_score \
        = _classify_crosslinks_by_score
