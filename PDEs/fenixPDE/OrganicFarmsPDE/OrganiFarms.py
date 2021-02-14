#We wanna solve u_t = D grad(u) +g*u -u^2, where u_x = 0 on x = 0, l, 1cm
from fenics import *
import matplotlib.pyplot as plt
T=2.0
num_steps = 10
dt = T/num_steps
D=1
tol = 1E-14
#define mesh and function space
mesh = IntervalMesh(50, 0, 1)
V = FunctionSpace(mesh, "CG", 2)

#BC
def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, Constant('1'), boundary)
#Initial Value
u_0 = Constant('1')
u_n = interpolate(u_0, V)
# Setting up the variational problem
u = TrialFunction(V)
v = TestFunction(V)
g = Expression('x[0] <= 0.5 + 1E-14 ? -1 : 1', degree=2)
F = u*v*dx+dt*D*dot(grad(u), grad(v))*dx - dt*v*g*u*dx + dt*u*u*v*dx -u_n*v*dx
a, L = lhs(F), rhs(F)
# solving the variational problem.
u = Function(V)
t=0
for n in range(num_steps):
    t += dt 
    
    solve( a == L, u, bc)

    u_n.assign(u)
# plotting solution

plot(u)
plt.show()
#Solve

#Data
