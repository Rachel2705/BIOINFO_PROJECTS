import matplotlib.pyplot as plt
import numpy as np

# Generate example data
x = np.linspace(0, 2 * np.pi, 100)
y_sin = np.sin(x)
y_cos = np.cos(x)
y_sum = y_sin + y_cos

# Create a plot with three subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

# Plot 1: Sine Wave
ax1.plot(x, y_sin, label='Sinusoidal', color='blue')
ax1.set_ylabel('Amplitude')
ax1.legend()

# Plot 2: Cosine Wave
ax2.plot(x, y_cos, label='Cosinusoidal', color='green')
ax2.set_ylabel('Amplitude')
ax2.legend()

# Plot 3: Sum of Sine and Cosine
ax3.plot(x, y_sin, label='Sinusoidal', color='blue', linestyle='dashed')
ax3.plot(x, y_cos, label='Cosinusoidal', color='green', linestyle='dashed')
ax3.plot(x, y_sum, label='Sum', color='red')
ax3.set_xlabel('X-axis')
ax3.set_ylabel('Amplitude')
ax3.legend()

# Title for the entire figure
fig.suptitle('Combined Sine and Cosine Waves')

# Adjust layout to prevent clipping of ylabel
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# Save the plot as an image
plt.savefig('interesting_plot.png')

# Show the plot
plt.show()
