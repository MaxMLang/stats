import numpy as np
import matplotlib.pyplot as plt

# Parameters for the Gaussian distributions
std_dev = np.sqrt(2)  # Standard deviation set to sqrt(2) to get variance 2
y_selected_values = [2, 4, 6]  # Selected y values
modes = [0, 2, 4]  # Modes for the Gaussian distributions at each y value (identical to the means for a Gaussian)

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Plot each Gaussian distribution vertically for the selected y values
for y, mode in zip(y_selected_values, modes):
    x = np.linspace(mode - 4 * std_dev, mode + 4 * std_dev, 300)
    y_gaussian = np.exp(-(x - mode)**2 / (2 * std_dev**2)) / (std_dev * np.sqrt(2 * np.pi))
    ax.plot(y_gaussian + y, x, color='blue')  # Shift the Gaussian to the y position

    # Add a point for the mean of the distribution
    ax.plot(y, mode, 'ro')  # 'ro' is the matplotlib code for a red circle

# Plot a linear dashed line through the red dots (means of the distributions)
ax.plot(y_selected_values, modes, 'k--')

# Set titles and labels
ax.set_title('Vertical Gaussian Distributions with Linear Trend through Means')
ax.set_xlabel('$y$')
ax.set_ylabel('Density')

# Show the plot
plt.show()

