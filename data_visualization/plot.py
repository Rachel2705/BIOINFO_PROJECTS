# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np

# generate example data
x = np.linspace(0, 10, 100)
y = np.sin(x)

# create a plot
plt.plot(x, y)
plt.title('Sinusoidal Function')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')

# save the plot as an image
plt.savefig('plot.png')
