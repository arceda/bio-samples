import numpy as np
from pylab import *

matrix = np.genfromtxt('c++/x60189_4_fitness_by_iteration.txt', delimiter=',')

x = matrix[:,0]
y = matrix[:,1]

plot(x, y)
xlabel('Iterations')
ylabel('Fitness')
title('Crow search in x60189_4')
draw()
show()