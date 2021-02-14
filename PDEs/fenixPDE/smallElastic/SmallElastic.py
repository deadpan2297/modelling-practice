#We want to solve the elasticity problem for a body Omega. This is
# -Div * sigma = f, in Omega
# sigma = lambda * tr(eps)I+2*mu*eps
# eps = (Div(u) + (Div(u))^T)/2
# For sigma the stress tensor, f the body force per unit volume, 
#lambda & mu  Lames parameters, I identity, tr trace, eps is the 
#strain-rate tensor, and u is just the displacement vector field.
#
#We will model a clamped beam deformed under its own weight in 3D. 
#This is f = (0,0,-rg), r the density and g gravity.  The beam has
# a square cross section with width W and Length L. u=u_D=(0,0,0)
#at the clamped end, x=0. The rest of the boundary is traction free.
#(T=0).

from fenics import *
from ufl import nabla_div
import matplotlib.pyplot as plt
#Variables Scaled
L=1; W=0.2
mu=1
rho=1
delta=W/L
gamma=0.4*delta**2
beta=1.25
lambda_=beta
g=gamma

#Mesh and functionspace
mesh = BoxMesh(Point(0,0,0), Point(L,W,W), 10, 3, 3)
V=VectorFunctionSpace(mesh, 'P', 1)

#BC
tol = 1E-14

def clamped_boundary(x,on_boundary):
    return on_boundary and x[0]<tol

bc = DirichletBC(V, Constant((0,0,0)), clamped_boundary)

#Strain and stress terms
def eps(u):
    return 0.5*(nabla_grad(u)+nabla_grad(u).T)

def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*eps(u)
#Variational PRoblem
u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Constant((0,0,-rho*g))
T = Constant((0,0,0))
a = inner(sigma(u),eps(v))*dx
L = dot(f,v)*dx + dot(T,v)*ds
#Solution
u = Function(V)
solve(a == L, u, bc)
#Plot
plot(u, title='Displacement', mode='displacement')

#Plotting stress
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)
von_Mises = sqrt(3./2*inner(s,s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)
plot(von_Mises, title='Stress intensity')

#Computer magnitude of displacement
u_magnitude = sqrt(dot(u,u))
u_magnitude = project(u_magnitude, V)
plot(u_magnitude, 'Displacement magnitude')

#vtkfile
File('smallElasticSolution/solution.pvd') << u
File('smallElasticSolution/von_mises.pvd') << von_Mises
File('smallElasticSolution/magnitude.pvd') << u_magnitude

plt.show()
