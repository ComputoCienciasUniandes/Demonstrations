import numpy as np
from numpy import sin, cos

from scipy.constants import g
import scipy.integrate as integrate

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from matplotlib.widgets import Slider

L = 5.0  # longitud del pendulo en m

th = 45.0 # th angulo inicial (grados)
w = 0.0 # velocidad angular inicial (grados/segundo) 

def analytic(th, w, t): 
    theta = th*cos(omega*t) + (w/omega)*sin(omega*t) # solucion analitica
    theta_dot = omega*(-th*sin(omega*t)+(w/omega)*cos(omega*t)) # velocidad (rad/s)
    return theta, theta_dot

def derivs(state, t):
    # solucion numerica
    dydx = np.zeros_like(state)
    dydx[0] = state[1]
    dydx[1] = -(omega**2)*sin(state[0])
    return dydx

dt = 0.1
t_max = 20
t = np.arange(0.0, t_max, dt) # arreglo de tiempos a evaluar
N = len(t)

# propiedades predeterminadas
a_color = 'grey'
alpha = 0.5

fig = plt.figure(figsize=(8, 4))

# cuadricula para las gráficas
gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0:, 0])
ax2 = plt.subplot(gs[0, 1:])
ax3 = plt.subplot(gs[1:, -1])

ax1.axis('equal')
ax1.set_xlim(-16, 16)
ax1.set_ylim(-16, 16)

axes = [ax1, ax2, ax3] # agrupación de los ejes
xlabels = ['$x$ (m)', 'Time (s)', r'$\theta (t)$ (rad/s)']
ylabels = ['$y$ (m)', r'$\theta (t)$ (rad/s)', r'$\dot{\theta} (t)$ (rad/s$^2$)']

# configura los ejes
for (i, ax) in enumerate(axes):
    ax.set_xlabel(xlabels[i])
    ax.set_ylabel(ylabels[i])
    ax.grid()

# lista de objetos de matplotlib, por eje, dinámicos y estáticos
ax1_dyna = [None]*2
ax1_dyna[1] = ax1.plot([], [], 'o-', lw = 2, color = a_color, alpha = alpha, label = r"$\sin\theta \approx \theta$")[0]
ax1_dyna[0] = ax1.plot([], [], 'o-', lw = 2, label = r"$\sin\theta$")[0]

ax2_static = ax2.plot([], [])[0], ax2.plot([], [], c = a_color, alpha = alpha)[0]
ax3_static = ax3.plot([], [])[0], ax3.plot([], [], c = a_color, alpha = alpha)[0]

ax2_dyna = ax2.plot([], [], 'o', c = 'b')[0], ax2.plot([], [], 'o', c = a_color, alpha = alpha)[0]
ax3_dyna = ax3.plot([], [], 'o', c = 'b')[0], ax3.plot([], [], 'o', c = a_color, alpha = alpha)[0]

fig.tight_layout()
fig.subplots_adjust(left=0.1, bottom=0.35)
ax1.legend(fontsize=12).get_frame().set_alpha(0.0) # legenda en la gráfica

# subejes para los sliders
axlength = plt.axes([0.25, 0.1, 0.65, 0.03]) 
axangle = plt.axes([0.25, 0.15, 0.65, 0.03])
axspeed = plt.axes([0.25, 0.20, 0.65, 0.03])

# sliders
slength = Slider(axlength, 'Length', 0.1, 10, valinit=L)
sangle = Slider(axangle, 'Angle', 0.1, 180.0, valinit=th)
sspeed = Slider(axspeed, 'Angular Speed', -50.0, 50.0, valinit=w)

def init(*val):
    """
    Función que inicializa la animación, calcula la solución a las ecuaciones.
    """
    global xs, ys, positions, velocities, omega
    
    L, th, w = slength.val, sangle.val, sspeed.val
    omega = np.sqrt(g/L)
    state = np.radians([th, w]) # estado inicial

    n_sol, n_speed = integrate.odeint(derivs, state, t).T # resuelve la ecuacion diferencial
    a_sol, a_speed = analytic(state[0], state[1], t) # solucion analitica
    
    positions, velocities = [n_sol, a_sol], [n_speed, a_speed]
    xlims = [(0, max(t)), (min(a_sol), max(a_sol))]
    ylims = [(min(a_sol), max(a_sol)), (min(a_speed), max(a_speed))]
    
    for i in range(2):
        axes[i+1].set_xlim(xlims[i])
        axes[i+1].set_ylim(ylims[i])
        ax2_static[i].set_data(t, positions[i])
        ax3_static[i].set_data(positions[i], velocities[i])
        
        ax1_dyna[i].set_data([], [])
        ax2_dyna[i].set_data([], [])
        ax3_dyna[i].set_data([], [])

    # cambio de coordenadas
    x1, y1 = L*sin(n_sol), -L*cos(n_sol)
    x2, y2 = L*sin(a_sol), -L*cos(a_sol)
    
    xs = [x1, x2]
    ys = [y1, y2]
    
    return tuple((ax1_dyna, ax2_dyna, ax3_dyna))

# redirecciona el evento
slength.on_changed(init)
sangle.on_changed(init)
sspeed.on_changed(init)

def animate(j):
    for i in range(2):
        ax1_dyna[i].set_data([0, xs[i][j]], [0, ys[i][j]])
        ax2_dyna[i].set_data(t[j], positions[i][j])
        ax3_dyna[i].set_data(positions[i][j], velocities[i][j])
    
    return tuple((ax1_dyna, ax2_dyna, ax3_dyna))

# fig.tight_layout()

ani = animation.FuncAnimation(fig, animate, np.arange(N),
                              interval=N/20.0, init_func=init)

# ani.save('Animation.gif', writer='imagemagick', fps = N/20.0, dpi = 50)
plt.show()
