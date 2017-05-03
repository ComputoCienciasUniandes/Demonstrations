import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

from astrohut import *
import time

def run_simulation(system, speeds, M, G, epsilon, tolerance):
    start = time.time()
    sim = Simulation(M, G, epsilon, tolerance = tolerance, pos = system, speeds = speeds)
    sim.start(0, 0.1, 0.01)
    return time.time() - start

def test_tolerance(N_tests = 10):
    N = 1096
    M_T = 5e10/1e7
    M = 2.0*M_T/N #two galaxies
    G = 44.97
    epsilon = 0.1

    N_first = int(N*0.5)

    system, speeds = example(N, M, G)
    tolerances = np.linspace(0, 1, N_tests)
    results = np.zeros((N_tests, 2))

    for i in range(N_tests):
        time_used = run_simulation(system, speeds, M, G, epsilon, tolerances[i])
        results[i] = tolerances[i], time_used

    return results

def test_particles(N_tests = 10):

    particles = np.logspace(1, 4, N_tests)
    results = np.zeros((N_tests, 3))

    for i in range(N_tests):
        N = int(np.ceil(particles[i]))
        if N%2 == 1:
            N += 1
        M_T = 5e10/1e7
        M = 2.0*M_T/N #two galaxies
        G = 44.97
        epsilon = 0.1

        system, speeds = example(N, M, G)
        time_used0 = run_simulation(system, speeds, M, G, epsilon, 0.0)
        time_used1 = run_simulation(system, speeds, M, G, epsilon, 1.0)
        results[i] = N, time_used0, time_used1

    return results

tolerance_results = test_tolerance()
particles_results = test_particles()
np.savetxt("tolerance_results.dat", tolerance_results)
np.savetxt("particles_results.dat", particles_results)

fig = plt.figure()
plt.plot(tolerance_results[:,0], tolerance_results[:,1], "-o")
plt.xlabel("Tolerance")
plt.ylabel("Time used (s)")
fig.savefig("Tolerance_results.pdf")

fig = plt.figure()
plt.plot(particles_results[:,0], particles_results[:,1], "-o", label="Tolerance = 0.0")
plt.plot(particles_results[:,0], particles_results[:,2], "-o", label="Tolerance = 1.0")
particles = np.logspace(1, 4)
expectedLOG = particles*np.log(particles)
expectedLOG *= particles_results[-1,2]/expectedLOG[-1]
expectedNOR = particles**2
expectedNOR *= particles_results[-1,1]/expectedNOR[-1]
plt.plot(particles, expectedLOG, label="$\mathcal{O}(N\logN)$")
plt.plot(particles, expectedNOR, label="$\mathcal{O}(N^2)$")
plt.xlabel("Number of particles")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Time used (s)")
plt.legend()
fig.savefig("Particles_results.pdf")
