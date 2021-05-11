from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

#for experiment E1
from BDMsolver_ExperimentE1 import *

#for experiment E2
#from BDMsolver_ExperimentE2 import *

def createPlotData(case):
    #load mesh and boundary data
    mesh = Mesh(case+".xml")
    boundaries = MeshFunction("size_t", mesh, case+"_facet_region.xml")

    #create solution of Stokes Problem
    z = Stokes(mesh, boundaries)
    w = z.solving()
    u,p= split(w)

    C = FunctionSpace(mesh, "DG", 0)

    #calculate shear stress tau
    M = as_matrix(((0.0, 1.0), (0.0, 0.0)))
    tau_xz = inner(2*z.mu(u)*sym(grad(u)), M)
    var_proj = project(tau_xz, C)
    
    #extract shear stress at the bottom
    daten = []
    dofmap = C.dofmap()     #map of dofs (dofmap is a class)
    for f in facets(mesh):      #for each facet of the mesh
        if boundaries[f.index()] == 1:      #check if the facet belongs to the bottom boundary
            c = list(cells(f))[0]       #list of cells at the facet
            c2d = dofmap.cell_dofs(c.index())       #acccessing dof indices of cell 
            v = var_proj.vector()[c2d][0]   #
            x = 0.0  #f.mid_point()
            for vert in vertices(f): # for each vertex at the facet (there are two vertices for each facet)
                x = x + vert.x(0)   # vertex_1 x-coordinate+ vertex_2 x-coordinate / 2
            x = x/2.0
            daten.append((x,v)) #append facet midpoint and the value
    daten.sort()
    return daten

##########create plot for 'oga' P2P1 implementation############
###########################################

#extract columns from file
def extract_col(file, col):
    column =[]
    with open(file) as f:
        for line in f:
            chars=[]
            line=line.split(' ')
            for char in line:
                if char not in ['',' ']:
                    chars.append(char)
            column.append(chars[col-1])
    return column

#data for Experiment E1
data = 'oga1e000.txt'

#data for Experiment E2
#data = 'oga1e001.txt'
xOga = extract_col(data,1)
u1Oga = extract_col(data,2)
u2Oga = extract_col(data,3)
tauOga = extract_col(data,4)
#calculate norm of u = (u1, u2)
vOga = []
for i in range(len(u1Oga)):
    vOga.append(sqrt(float(u1Oga[i])**2 + float(u2Oga[i])**2))

#plot tau for BDM solver an different mesh resolutions 
##    ATTENTION  mesh size 5 is computationally expensive
tau_data20 = createPlotData('mesh_E_20')
plt.plot([x[0]/1000.0 for x in tau_data20],[x[1]/1000.0 for x in tau_data20], label = 'BDM; h=20m')
tau_data10 = createPlotData('mesh_E_10')
plt.plot([x[0]/1000.0 for x in tau_data10],[x[1]/1000.0 for x in tau_data10], label = 'BDM; h=10m')
#tau_data5 = createPlotData('mesh_E_5')
#plt.plot([x[0]/1000.0 for x in tau_data5],[x[1]/1000.0 for x in tau_data5], label = 'BDM; h=5m')

#plot 'oga' data
plt.plot([float(x)*5 for x in xOga], [float(tau) for tau in tauOga], ':', label = 'P2P1-Stokes')
plt.xlabel("x in km")
plt.ylabel(r"$\tau_{xz}(z_b)$ in (kPa)")
plt.legend()
plt.show()
