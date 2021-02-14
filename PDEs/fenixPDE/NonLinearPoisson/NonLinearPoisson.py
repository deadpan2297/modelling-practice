#We are solving -Div(q(u)Div(u))=f

from fenics import *

def q(u):
    return 1+ u**2


#make the mesh and function space
mesh = UnitSquareMesh(8,8)
V = FunctionSpace(mesh, 'P', 1)
u_D = Expression('x[0]+2*x[1]+1', degree=1)

#BC
def boundary(x, on_boundary):
    return on_boundary

bc=DirichletBC(V,u_D, boundary)

#Vtkfile
vtkfile = File('nonlinearPoison/solution.pvd')

#Solving
u=Function(V)
v=TestFunction(V)
f=Expression('-10*x[0]-20*x[1]-10', degree=1)
#The variational problem
F=q(u)*dot(grad(u),grad(v))*dx - f*v*dx 
solve(F==0, u, bc)
vtkfile<<u
#Output
plot(u)


