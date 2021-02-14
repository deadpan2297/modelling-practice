import matplotlib.pyplot as plt
from fenics import *
import numpy as np

T=100.0
num_steps = 50
dt = T/num_steps
alpha=3
beta=1.2

#Mesh and function space
nx=ny=8
mesh  = UnitSquareMesh(nx,ny)
V=FunctionSpace(mesh,'P',1)

#BC

def boundary(x, on_boundary):
    return on_boundary

u_D = Expression('1+x[0]*x[0] + alpha*x[1]*x[1] + beta*t',degree=2,alpha=alpha, beta=beta, t=0)

bc = DirichletBC(V, u_D, boundary)

#initial value

u_n = interpolate(u_D, V)

#Variational problem
u = TrialFunction(V)
v = TestFunction(V)
f=Constant(beta - 2 - 2*alpha)

F=u*v*dx + dt*dot(grad(u),grad(v))*dx-(u_n+dt*f)*v*dx
a,L=lhs(F),rhs(F)
#Saving solition to VTK
vtkfile=File('heat_gaussian/solution.pvd')

#Time stepping
u=Function(V)
t=0
for n in range(num_steps):

    t +=dt
    u_D.t = t
    solve(a==L, u ,bc)
    vtkfile<<(u,t)

    plot(u)

    u_n.assign(u)
plt.show()
