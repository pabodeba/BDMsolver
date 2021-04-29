from dolfin import *

class Stokes():
    #set relevant physical quantities
    rho = Constant(910.0)   #ice density in kg m^-3 
    g = as_vector([0.,-9.81])   #graviational acceleration 2D
    beta = Constant(1e8) #basal slipperiness
    lmbda = Constant(10.0) #penaltyterm


    def __init__(self, mesh, boundaries):
        self.mesh = mesh
        self.boundaries = boundaries
        self.V = FiniteElement("BDM", self.mesh.ufl_cell(), 1)
        self.Q = FiniteElement("DG", self.mesh.ufl_cell(), 0)

        self.VQ = FunctionSpace(self.mesh, self.V*self.Q)
        #impenetrability boundary condition as a dirichlet condition
        self.dirichletbc= DirichletBC(self.VQ.sub(0), Constant((0.0, 0.0)), self.boundaries, 1)

        self.dx = Measure('dx', domain=self.mesh)
        self.ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        self.h = CellDiameter(self.mesh) #max. diameter of the elements
        self.n = FacetNormal(self.mesh) # facet normal

    def mu(self,u):  #effective viscosity
        A_rate = Constant(pow(10,-16))   #ice flow rate in SI units
        eps_e = sqrt(0.5*inner(sym(grad(u)), sym(grad(u)))+1e-4)
        return 0.5* pow(A_rate,-1.0/3.0)*pow(eps_e, 1.0/3.0-1.0) # 3.0 is the exponent d in Glenn's flow law
    
    def sigma(self, u, p):    #Stress tensor
        return 2.0*sym(grad(u)) - p*Identity(u.geometric_dimension())

    #bilinear forms for the discontinuous weak formulation interior penalty method)
    def A(self,u, v):
        su = self.mu(u)*2.0*sym(grad(u))
        sv = self.mu(u)*2.0*sym(grad(v))
        tria = inner(su, grad(v))*self.dx
        edge = -inner(avg(su), 2*avg(outer(self.n, v)))*dS
        symm = -inner(avg(sv), 2*avg(outer(self.n, u)))*dS
        penal = Stokes.lmbda/avg(self.h) * inner(avg(self.mu(u))*2*avg(outer(self.n, u)),2*avg(outer(self.n, v)))*dS
        return tria + edge + symm + penal

    def B(self, v,p):
        tria = -inner(p, div(v))*self.dx 
        edge =  + inner(avg(p),2*avg(dot(v,self.n)))*dS   
        return tria + edge

    def C(self,u,v):
        return -Stokes.beta*inner(u,v)*self.ds(1)    

    def L(self,v):
        return Stokes.rho*inner(Stokes.g,v)*self.dx 
    
    def solving(self):
        w1 = Function(self.VQ)
        u, p = split(w1)
        w2 = TestFunction(self.VQ) 
        v, q = split(w2)

        F = self.A(u,v) + self.B(v,p) + self.B(u,q) +self.C(u,v) - self.L(v)

        solve(F == 0, w1, self.dirichletbc)
        print("Done solving")
        return w1
