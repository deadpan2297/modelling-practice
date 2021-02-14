from ddeint import ddeint
from numpy import *
import matplotlib.pyplot as plt

samples = 1000
def model(Y,t):
    return -Y(t-3*cos(Y(t))**2)

def history(t):
    return 1


ts = linspace(0,150, samples)
ys = ddeint(model, history, ts)
ts1 = linspace(0,150, 5*samples)
ys1 = ddeint(model, history, ts1)
ts2 = linspace(0,150, 10*samples)
ys2 = ddeint(model, history, ts2)

fig, ax = plt.subplots()

plt.xlabel("y(t)")
plt.ylabel("t")
plt.title('Simple DDE with varying sample size')

l1,=ax.plot(ts,ys)
l2,=ax.plot(ts1,ys1)
l3,=ax.plot(ts2,ys2)

ax.legend((l1,l2,l3), ('%s Samples'%(samples), '%s Samples'%(5*samples), '%s Samples'%(10*samples)), loc='upper right', shadow=True)
ax.set_xlabel("t")
ax.set_ylabel("y(t)")
plt.show()
