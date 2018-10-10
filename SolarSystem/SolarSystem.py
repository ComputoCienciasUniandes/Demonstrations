import numpy as np
from numpy import sin, cos

from datetime import datetime

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

ref_time = datetime(1970,1,1, 0, 0)
init_time = datetime(2016, 12, 11, 0, 0)
init_time = (init_time - ref_time).total_seconds()

year_to_seconds = 365.25*24*60*60

data = np.genfromtxt("coordinates.csv", delimiter=",")

masses = data[:, 1]/data[0, 1]
positions = data[:, 2:5]
speeds = data[:, 5:]

G = 4*(np.pi**2) 

def acceleration(R):
    a = np.zeros((10, 3))
    for i in range(10):
        for j in range(10):
            d = R[j] - R[i]
            if i != j:
                mag = np.sqrt(d.dot(d))
                mag = mag**3
                
                a[i] += G*masses[j]*d/mag    
    return a

def solver(positions, speeds, t, dt):
    N = int(t/dt)
    positions_in_time = np.zeros((N, 10, 3))
    x = positions.copy()
    v = speeds.copy()
    positions_in_time[0] = x
    
    for i in range(N-1):
        v_half = v + 0.5*dt*acceleration(x)
        x += dt*v_half
        v = v_half + 0.5*dt*acceleration(x)   
        positions_in_time[i+1] = x
        
    return positions_in_time

T_max = 252
dt = 10/365.25
positions_in_time = solver(positions, speeds, T_max, dt)

labels = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto :("]
colors = ["yellow", "grey", "orange", "blue", "red", "orange", "yellow", "cyan", "blue", "black"]

step = 50
shorten = positions_in_time[::step]
N = len(shorten)

fig = plt.figure()
ax = p3.Axes3D(fig)

text = ax.text(-40, -30, 5, "")
fixed = [ax.plot(shorten[:, i, 0], shorten[:, i, 1], shorten[:, i, 2], c = colors[i], lw = 1) for i in range(10)]
plots = [ax.plot([], [], [], "o", label = labels[i], color = colors[i])[0] for i in range(10)]
plots[0].set_marker("*")
plots[0].set_markersize(20)

ax.set_xlabel("$x$ (AU)")
ax.set_ylabel("$y$ (AU)")
ax.set_zlabel("$z$ (AU)")

plt.legend(numpoints=1)

def init():
    for (j, line) in enumerate(plots):
        line.set_data([], [])
        line.set_3d_properties([])
    text.set_text("")
    return plots, text

def update(i):
    for (j, line) in enumerate(plots):
        line.set_data(shorten[i, j, 0], shorten[i, j, 1])
        line.set_3d_properties(shorten[i, j, 2])
    time = i*dt*year_to_seconds*len(positions_in_time)/float(N) + init_time
    time = datetime.utcfromtimestamp(time)
    text.set_text(time.strftime('%Y-%m-%d'))
    return plots, text

ani = animation.FuncAnimation(fig, update, N)
#ani.save("Planets.gif", writer = "imagemagick", fps = N/30, dpi = 50)
plt.show()
