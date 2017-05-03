from astrohut import *

N = 4096
M_T = 5e10/1e7
M = 2.0*M_T/N #two galaxies
G = 44.97
epsilon = 0.01

galaxy1, galaxy2, system, speeds = example(N, M, G)

sim = Simulation(M, G, system, speeds, epsilon = epsilon,
        tolerance = 1.0, threads = -1)
sim.start(0, 10.0, 0.01)

data = read_output()
ani = animate(data, N)
plt.show()
#ani.save("approx.mp4", writer = "ffmpeg", fps = 24, dpi = 120)
