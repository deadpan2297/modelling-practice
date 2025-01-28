import matplotlib.pyplot as plt
from fenics import *
import numpy as np

T=2
num_steps = 250
dt = T/num_steps
alpha=3
beta=1.2

#Mesh and function space
nx=ny=30
mesh  = RectangleMesh(Point(-2,-2), Point(2,2), nx, ny)
V=FunctionSpace(mesh,'P',1)

#BC
def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, Constant(0), boundary)

#initial value
u_0 = Expression('exp(-a*pow(x[0],2) - a*pow(x[1],2))', degree=2, a=5)
u_n = interpolate(u_0, V)

#Variational problem
u = TrialFunction(V)
v = TestFunction(V)
f=Constant(0)

F=u*v*dx + dt*dot(grad(u),grad(v))*dx-(u_n+dt*f)*v*dx
a,L=lhs(F),rhs(F)
#Saving solition to VTK
vtkfile=File('output2/GaussianSolution.pvd')

#Time stepping
u=Function(V)
t=0
for n in range(num_steps):

    #Update current time
    t +=dt
    #compute the solution
    solve(a==L, u ,bc)
    #Save to file and plot
    vtkfile<<(u,t)
    plot(u)
    #Update previous solutio:
    u_n.assign(u)
plt.show()
