#This is an script using ddeint to show the
#hopf bifurcation of the DDE
# y' = -y(t-tau)
# at the first point bifurcation = pi/2.
# Note that there is a Hopf bifurcation at
# tau = (-1)^n (pi/2 + n*pi), n = ..., -2, -1, 0, 1, 2, ...
#
import matplotlib.pyplot as plt
from ddeint import ddeint
import numpy as np
#Delays
tau1 = np.pi/2+0.1
tau2 = np.pi/2
tau3 = np.pi/2-0.1
#Model
def model(Y,t,tau):
    return -Y(t-tau)
#History Function
def history(t):
    return 1
#Simulate the DDE with 3 different delays
ts = np.linspace(0,100,5000)
ys1 = ddeint(model, history, ts, fargs=(tau1,))
ys2 = ddeint(model, history, ts, fargs=(tau2,))
ys3 = ddeint(model, history, ts, fargs=(tau3,))

#Plotting 3 Subplots in 1 major plot
fig, ax = plt.subplots(3,1,constrained_layout=True)
#Set x-axis and title
plt.xlabel('t')
fig.suptitle(r'Solution to  $\frac{dy}{dt} = -y(t-\tau)$')
#Plot each simulation on the subplots
l1 = ax[0].plot(ts,ys1, color='red')
l2 = ax[1].plot(ts,ys2, color='orange')
l3 = ax[2].plot(ts,ys3)
#Set title of each subplot
ax[0].set_title(r'$\tau = \frac{\pi}{2} + 0.1$')
ax[1].set_title(r'$\tau = \frac{\pi}{2}$')
ax[2].set_title(r'$\tau = \frac{\pi}{2} - 0.1$')
#Display the plot
plt.show()
