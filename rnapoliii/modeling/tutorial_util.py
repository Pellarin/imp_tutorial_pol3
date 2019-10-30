import IMP.pmi.io.crosslink
import IMP.pmi.macros
import IMP.pmi.dof
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


# Monkey patch DegreesOfFreedom._setup_srb to incorporate a post-2.11 bug fix
def _setup_srb(self,hiers,max_trans,max_rot,axis):
    if axis is None:
        srbm = IMP.pmi.TransformMover(hiers[0][0].get_model(), max_trans, max_rot)
    else:
        srbm = IMP.pmi.TransformMover(hiers[0][0].get_model(),axis[0],axis[1],max_trans, max_rot)
    srbm.set_was_used(True)
    super_rigid_rbs,super_rigid_xyzs = IMP.pmi.tools.get_rbs_and_beads(hiers)
    ct = 0
    self.movers_particles_map[srbm]=[]
    for h in hiers:
        self.movers_particles_map[srbm]+=IMP.atom.get_leaves(h)
    for xyz in super_rigid_xyzs:
        srbm.add_xyz_particle(xyz)
        ct+=1
    for rb in super_rigid_rbs:
        srbm.add_rigid_body_particle(rb)
        ct+=1
    if ct>1:
        return srbm
    else:
        return 0



if IMP.__version__ == '2.11.1':
    IMP.pmi.macros.AnalysisReplicaExchange.align = _align
    IMP.pmi.io.crosslink.CrossLinkDataBase.append_database = _append_database
    IMP.pmi.io.crosslink.CrossLinkDataBase.classify_crosslinks_by_score \
        = _classify_crosslinks_by_score
    IMP.pmi.dof.DegreesOfFreedom._setup_srb = _setup_srb
