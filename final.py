import pychrono as chrono
import pychrono.fea as fea
import pychrono.pardisomkl as mkl
import pychrono.irrlicht as chronoirr
import errno
import os
import copy


# The path to the Chrono data directory containing various assets (meshes, textures, data files)
# is automatically set, relative to the default location of this demo.
# If running from a different directory, you must change the path to the data directory with: 
#chrono.SetChronoDataPath('path/to/data')


#Material properties for material 1 (Stem)
E1 = 16e9 # rubber 0.01e9, steel 200e9
v1 = 0.3
dampingCoeff1 = 0.1
density1 = 400

#Material properties for material 2 (Leaf)
E2 = 0.01e9 # rubber 0.01e9, steel 200e9
v2 = 0.3
dampingCoeff2 = 0.1 #make it higher (*10)
density2 = 10

timestep = 1.0
finalTime = 50

writeVTU = True
writeTXT = True

# Create a Chrono::Engine physical system
my_system = chrono.ChSystemSMC()

# Create a mesh
mesh = fea.ChMesh()
my_system.Add(mesh)

## Create a mesh, that is a container for groups
## of elements and their referenced nodes.

# Create some nodes (with default mass 0)
nodeA = fea.ChNodeFEAxyzrot(chrono.ChFrameD(chrono.ChVectorD(0, 0, 0)))
nodeB = fea.ChNodeFEAxyzrot(chrono.ChFrameD(chrono.ChVectorD(1, 0, 0)))
nodeC = fea.ChNodeFEAxyzrot(chrono.ChFrameD(chrono.ChVectorD(0, 1, 0)))
nodeD = fea.ChNodeFEAxyzrot(chrono.ChFrameD(chrono.ChVectorD(0, 0, 1)))

nodeA.SetYoungModulus(E1)
nodeB.SetYoungModulus(E2)
nodeA.SetDensity(density1)
nodeB.SetDensity(density2)
nodeA.SetPoissonRatio(v1)
nodeB.SetPoissonRatio(v2)
nodeA.SetBeamRaleyghDamping(dampingCoeff1)
nodeB.SetBeamRaleyghDamping(dampingCoeff2)

mesh.AddNode(nodeA)
mesh.AddNode(nodeB)
mesh.AddNode(nodeC)
mesh.AddNode(nodeD)

element1= fea.ChElementTetrahedron()
element2= fea.ChElementTetrahedron()


element1.SetNodes(nodeA, nodeB, nodeC, nodeD)
#element2.SetNodes(nodeE, nodeF, nodeG, nodeH)

mesh.AddElement(element1)

# We want gravity effect on FEA elements in this demo
mesh.SetAutomaticGravity(True)
my_system.Set_G_acc(chrono.ChVectorD(0,-9.81, 0))



