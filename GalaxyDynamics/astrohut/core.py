import os
import ctypes
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

class Simulation:
    def __init__(self, M, G, pos, speeds, epsilon = 0.01, tolerance = 1,
                        output="Data/", print_energy = False, threads = -1):
        self.M = M
        self.G = G
        self.epsilon = epsilon
        self.tolerance = tolerance
        self.pos = pos
        self.speeds = speeds
        self.N = self.pos.shape[1]
        self.threads = threads
        self.print_energy = int(print_energy)
        """
        c init
        """
        self.c_double_p = ctypes.POINTER(ctypes.c_double)
        self.x = self.pos[0].ctypes.data_as(self.c_double_p)
        self.y = self.pos[1].ctypes.data_as(self.c_double_p)
        self.z = self.pos[2].ctypes.data_as(self.c_double_p)
        self.vx = self.speeds[0].ctypes.data_as(self.c_double_p)
        self.vy = self.speeds[1].ctypes.data_as(self.c_double_p)
        self.vz = self.speeds[2].ctypes.data_as(self.c_double_p)

        cwd = os.path.abspath(__file__)[:-len('core.py')]
        self.path = glob("%sbruteforce*.so"%cwd)[0]

        self.init_lib()
        self.output = output

        if os.path.isdir(self.output):
            files = glob("%s%s"%(self.output, "*.dat"))
            list(map(os.remove, files))
        else:
            os.makedirs(self.output)

    def init_lib(self):
        self.lib = ctypes.CDLL(self.path)

        self.lib.init_conditions.argtypes = (ctypes.c_int, ctypes.c_double,
                            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int)
        self.lib.solver.argtypes = (ctypes.c_double, ctypes.c_double,
                            ctypes.c_double, ctypes.c_char_p)

        self.lib.init_conditions(self.N, self.M, self.G, self.epsilon,
                            self.tolerance, self.print_energy, self.threads)
        self.lib.init_from_ram.argtypes = (self.c_double_p, self.c_double_p,
                self.c_double_p, self.c_double_p, self.c_double_p, self.c_double_p)

        self.lib.init_from_ram(self.x, self.y, self.z, self.vx, self.vy, self.vz)

    def start(self, t0, tmax, dt):
        if self.lib == None:
            self.init_lib()
        self.lib.solver(t0, tmax, dt, self.output.encode())
        del self.lib
        self.lib = None

    def box(self):
        self.lib.init_tree()

class Galaxy(object):
    def __init__(self, N_particles, diameter):
        self.N_particles = N_particles
        self.diameter = diameter
        self.positions = None
        self.center_of_mass = None
        self.distances = None
        self.create_positions()
        self.create_distances()

    def create_positions(self):
        self.positions = np.zeros((3, self.N_particles))
        theta = 2*np.pi*np.random.random(self.N_particles)
        r = 0.5*self.diameter*np.random.normal(size = self.N_particles)
        self.positions[0] = r*np.cos(theta)
        self.positions[1] = r*np.sin(theta)
        self.positions[2] = 0.02*self.diameter*2.0*(np.random.random(self.N_particles)-0.5)

        self.center_of_mass = np.mean(self.positions, axis = 1)
        self.positions = np.array([self.positions[i] - self.center_of_mass[i]\
                                                            for i in range(3)])

    def create_distances(self):
        self.distances = np.sum(self.positions**2, axis = 0)**0.5

def read_output(output="Data/"):
    files = glob("%s*_instant.dat"%output)
    for file_name in files:
        with open(file_name, "r") as reading:
            text = reading.read()
        text = text.replace(',', '.')
        with open(file_name, "w") as writing:
            writing.write(text)
    return [np.genfromtxt("%s%d_instant.dat"%(output, i)) for i in range(len(files))]

def speeds_generator(galaxy, G, m = 1):
    distances = galaxy.positions
    r = galaxy.distances
    theta = np.arctan2(distances[1, :], distances[0, :])
    n = len(r)
    s = np.zeros((3, n))
    for i in range(n):
        my_pos = distances[:, i]
        others_pos = np.array([(distances[0, :] - my_pos[i])**2 for i in range(3)]).sum(axis=0)
        # speed_mag = np.sqrt(G*m*(1/np.sqrt(others_pos)).sum())
        pos = np.where(r < r[i])[0]
        M = m*pos.shape[0]
        theta = np.arctan2(distances[1, i], distances[0, i])
        if M != 0:
            mag = np.sqrt(G*M/r[i])
            s[:2, i] = -r[i]*np.sin(theta), r[i]*np.cos(theta)
            s[:2, i] *= 1/np.sqrt(np.dot(s[:, i], s[:, i]))
            s[:2, i] *= mag
            # if speed_mag > mag:
            #     s[2, i] = (speed_mag - mag)*(np.random.random() - 0.5)
        else:
            s[:, i] = 0
    return 0.93*s

def rotation_matrixes(theta):
    x = np.array([[1, 0, 0], [0, np.cos(theta), - np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
    y = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
    z = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    return x, y, z

def example(N, M, G):
    N_first = int(0.5*N)
    N_second = int(N-N_first)
    theta = np.pi/4

    rotx, roty, rotz = rotation_matrixes(theta)

    galaxy1 = Galaxy(N_first, 55)
    galaxy2 = Galaxy(N_second, 55)

    speeds1 = speeds_generator(galaxy1, G, M)
    speeds2 = speeds_generator(galaxy2, G, M)

    # galaxy2.positions = np.dot(rotx, galaxy2.positions)
    # speeds2 = np.dot(rotx, speeds2)

    displacement = 2.5*55
    system = np.zeros((3, N))
    system[:, :N_first] = galaxy1.positions
    system[:, N_first:] = galaxy2.positions + displacement
    speeds = np.zeros_like(system)

    speeds[:, :N_first] = speeds1
    speeds[:, N_first:] = speeds2

    return galaxy1, galaxy2, system, speeds

def systemEnergy(directory = "Data/"):
    files = glob("%s*_energy.dat"%directory)
    energies = [np.genfromtxt("%s%d_energy.dat"%(directory, i)) for i in range(len(files))]
    total = np.array([energy.sum() for energy in energies])
    total = (total - total[0])/total[0]
    return total


"""
plotting
"""
def animate(data, N, dt = 0.01):
    fig = plt.figure(figsize=(16,9))
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.set_aspect("equal")
    plot1 = ax1.plot([], [], [], "o", ms=0.5, c="g", alpha = 0.5)[0]
    plot2 = ax1.plot([],[],[], "o", ms=0.5, c="b", alpha = 0.5)[0]
    ax1.set_xlabel("$x$ (kpc)")
    ax1.set_ylabel("$y$ (kpc)")
    ax1.set_zlabel("$z$ (kpc)")
    ax2 = fig.add_subplot(1, 2, 2)
    plot2d1 = ax2.plot([], [], "o", ms=0.5, c="g", alpha = 0.5)[0]
    plot2d2 = ax2.plot([],[], "o", ms=0.5, c="b", alpha = 0.5)[0]
    ax2.set_xlabel("$y$ (kpc)")
    ax2.set_ylabel("$z$ (kpc)")
    fig.suptitle('Time: 0')
    N_first = int(0.5*N)

    def update(i):
        temp = data[i]
        plot1.set_data(temp[:N_first,0], temp[:N_first,1])
        plot1.set_3d_properties(temp[:N_first, 2])
        plot2.set_data(temp[N_first:,0], temp[N_first:,1])
        plot2.set_3d_properties(temp[N_first:, 2])

        plot2d1.set_data(temp[:N_first,1], temp[:N_first,2])
        plot2d2.set_data(temp[N_first:,1], temp[N_first:,2])
        fig.suptitle('Time: %.2f Gy'%(i*dt))

    min_value, max_value = -100, 200
    ax1.set_xlim(min_value, max_value)
    ax1.set_ylim(min_value, max_value)
    ax1.set_zlim(min_value, max_value)
    ax2.set_xlim(min_value, max_value)
    ax2.set_ylim(min_value, max_value)
    ani = FuncAnimation(fig, update, frames=len(data), interval=0.01)
    return ani
