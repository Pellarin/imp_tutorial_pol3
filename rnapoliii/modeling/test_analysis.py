from __future__ import print_function

import IMP
import IMP.core
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.tools

import IMP.pmi.macros
import IMP.pmi.topology

# Hot fixes correcting minor bugs in IMP 2.11.1
import tutorial_util

import os
import sys

# Initialize IMP model
m = IMP.Model()

nmodels=10
rmffile="./results/150_xl_cryoem_3.rmf"

import IMP.pmi.output

hh=IMP.pmi.output.StatHierarchyHandler(m,rmffile)

#Total number of frames
print("Frames",len(hh))

# Describe the content of the first frame of the rmf file
print(hh[0])

#list down all the features names
for k in hh[0].features.keys(): print(k)
    

# For instance we can compute the distance between two residues

%pylab inline

p0=IMP.atom.Selection(hh,molecule="C31",residue_index=10).get_selected_particles()[0]
p1=IMP.atom.Selection(hh,molecule="C34",residue_index=10).get_selected_particles()[0]

d0=IMP.core.XYZ(p0)
d1=IMP.core.XYZ(p1)

#note that hh can be used as a list
plot([IMP.core.get_distance(d0,d1) for h in hh]);

figure()

# Or we can get the radius of gyration of the whole complex
ps=IMP.atom.Selection(hh).get_selected_particles()
plot([IMP.atom.get_radius_of_gyration(ps) for h in hh])

# To reduce I/O, we can store the data structure internal to hh, 
# so that it is not read directly from the files
# and it is faster

data=hh.data

# Then we plot the scores
plot([x.score for x in data])

figure() 

# finally we plot the distance of two crosslinked residues
plot([float(x.features["CrossLinkingMassSpectrometryRestraint_Score_|XL|29.APO.1|C31|91|C160|1458|0|CLASS_0|"]) for x in data]);


scores={}

for xl in xl1: 
    if not xl['IntraRigidBody']:
        scores.update({xl['XLUniqueSubID']:float(xl['IDScore'])})
    

x=[]
y=[]
for k in data[0].features.keys():
    if "Distance" in k:
        id=k.split("|")[2]
        if id in scores:
            x.append(scores[id])
            y.append(float(data[2].features[k]))

scatter(x,y);


are=IMP.pmi.macros.AnalysisReplicaExchange(m,
                 [rmffile],
                 best_models=nmodels,
                 alignment=False)

print(are)


are.set_rmsd_selection(molecules=["C31","C34","C53","C37","C82"])

are.cluster(10.0)

# see the contant of the "are" object
print(are)

#print the cluster info
for cluster in are:
    print(cluster)

for member in are[0]:
    print(member)
    
are.save_coordinates(are[0])

# slow!
are.plot_rmsd_matrix("rmsd_matrix.pdf")


keys=[k for k in are[0][0].features.keys() if "CrossLinkingMassSpectrometryRestraint" in k]
figure()
for ind in ["0","1","2"]:
    distances={}
    psi=[]
    sigma=[]
    for member in are[0]:
        for k in keys:
            if "CLASS_"+str(ind) in k and "Distance" in k:
                if k in distances:
                    distances[k].append(float(member.features[k]))
                else:
                    distances[k]=[float(member.features[k])]
            if "CLASS_"+str(ind) in k and "Psi" in k and "MonteCarlo" not in k:
                psi.append(float(member.features[k]))
            if "SIGMA" in k and "Sigma" in k and "MonteCarlo" not in k:
                sigma.append(float(member.features[k]))

    x=[distances[k] for k in distances]
    pylab.boxplot(x);
    figure()
    pylab.plot(range(len(psi)),psi);
    figure()
pylab.plot(range(len(sigma)),sigma);




### build the seed 

are.write_seed("seed.rmf3", 48)

are.compute_cluster_center(cluster=are[0])
print(are.precision(cluster=are[0]))
print(are.bipartite_precision(cluster1=are[0],cluster2=are[1]))

rmsf1=are.rmsf(cluster=are[0],molecule='C31');
plot(list(rmsf1.keys()),list(rmsf1.values()),marker=".",linewidth=0)
figure()

rmsf2=are.rmsf(cluster=are[0],molecule='C34');
plot(list(rmsf2.keys()),list(rmsf2.values()),marker=".",linewidth=0)

for mol in ['ABC23','ABC10beta','ABC14_5','ABC27','C25','AC40','C160','ABC10alpha','C128','AC19','C11','C17','C31','C34','C53','C37','C82']: 
    are.rmsf(cluster=are[0],molecule=mol);
ch1=IMP.pmi.tools.ColorHierarchy(are.stath1)
ch1.color_by_uncertainty()
are.save_coordinates(are[0])

density_names={'core': ['ABC23','ABC10beta','ABC14_5','ABC27','C25','AC40','C160','ABC10alpha','C128','AC19','C11','C17'],
               'C53': ['C53'], 
               'C37': ['C37'], 
               'C34': ['C34'], 
               'C82': ['C82'], 
               'C31': ['C31']}

# you can iterate on the clusters
for n,a in enumerate(are):
    are.save_densities(cluster=a,density_custom_ranges=density_names,prefix="Cluster-"+str(n))

# it is slow
are.contact_map(cluster=are[0]);









