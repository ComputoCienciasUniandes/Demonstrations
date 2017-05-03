from astrohut import *

N = 4096
data = read_output()
ani = animate(data, N)
plt.show()
