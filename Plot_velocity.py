from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

#for experiment E1
from BDMsolver_ExperimentE1 import *

#for experiment E2
#from BDMsolver_ExperimentE2 import *

################Create Experiment solution#################
############################################################
#Load meshe and boundary
case = "ExperimentE_20"
GeoFolder = "./Geometries"

mesh = Mesh(GeoFolder+"/"+case+".xml")
boundaries = MeshFunction("size_t", mesh, GeoFolder+"/"+case+"_facet_region.xml")

z = Stokes(mesh, boundaries)
w = z.solving()
u,p = split(w)
#V = z.VQ.sub(0).collapse()
V = VectorFunctionSpace(mesh, 'CG', 2)
u_proj = project(u, V)

#####################Create data for plots##########################
###########################################################

def PlotFunctionOfX(u_cg, V):
    dofs = []
    cell_to_facets = mesh.topology()(2,1) #get mesh topology, returns the topology object associated with the mesh.
    for cell in cells(mesh):        #for each cell in the mesh
        facets = cell_to_facets(cell.index()) #for the edges of the cell
        for facet in facets:
            if boundaries[facet] == 2:      #if its a surface edge
                #Return the dof indices associated with the closure of facets and facet indices
                dofs_ = V.dofmap().entity_closure_dofs(mesh, 1, [facet]) 
                for dof in dofs_:
                    dofs.append(dof)
    dofs.sort()                

    unique_dofs = np.array(list(set(dofs)), dtype=np.int32)
    boundary_coords = V.tabulate_dof_coordinates()[unique_dofs]
    data = []

    for i, dof in enumerate(unique_dofs):
        temp_list = [boundary_coords[i][0], u_cg.vector()[dof]]
        data.append(temp_list)
    data.sort()
    
    xValues = []
    uValues = []
    i=0
    #calculate norm of velocity vector
    while 2*i <len(data):
        xValues.append(data[2*i][0])
        uValues.append(sqrt(data[2*i][1]**2 + data[2*i+1][1]**2))
        i=i+1
    
    return xValues, uValues

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
#calculate norm of u = (u1, u2)
vOga = []
for i in range(len(u1Oga)):
    vOga.append(sqrt(float(u1Oga[i])**2 + float(u2Oga[i])**2))


xValues, uValues = PlotFunctionOfX(u_proj, V)
plt.plot([x/1000.0 for x in xValues], uValues, label = 'BDM-Stokes')

plt.plot([float(x)*5 for x in xOga], vOga, ':', label = 'P2P1-Stokes')
plt.xlabel("x in km")
plt.ylabel("u ($\mathregular{m a^{-1}}$)")
plt.legend()
plt.show()