##This is just an example program from the book (Solving PDEs in Python)
from fenics import * 
import matplotlib.pyplot as plt

#Create the mesh and define function space
mesh = UnitSquareMesh(8,8)
V= FunctionSpace(mesh, 'P', 1)

#Boundary conditions
u_D = Expression('1-x[0]*x[0]+2*x[1]*x[1]', degree = 2)

def boundary(x, on_boundary):
    return on_boundary


bc=DirichletBC(V,u_D, boundary)

#The variational problem a(u,v)=L(v)
u=TrialFunction(V)
v=TestFunction(V)
f=Constant(-6.0)
a=dot(grad(u),grad(v))*dx
L=f*v*dx

#Computing the solution
u=Function(V)
solve(a==L, u, bc)

#Plotting
plot(u)
plot(mesh)

#Save the solution to file in VTK format
vtkfile = File('solution.pvd')
vtkfile << u

plt.show()
