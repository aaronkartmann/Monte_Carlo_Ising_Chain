import numpy as np


def bin_data(data, bin_size):
    """Split data into consecutive bins and return the mean of each bin."""
    data = np.asarray(data)

    if data.ndim != 1:
        raise ValueError("data must be one-dimensional.")
    if bin_size < 1:
        raise ValueError("bin_size must be positive.")

    num_bins = len(data) // bin_size # Calculate the number of complete bins that can be formed with the given bin size

    if num_bins < 2:
        raise ValueError("At least two complete bins are required.")

    data = data[:num_bins * bin_size] # Trim data to fit exactly into complete bins
    return data.reshape(num_bins, bin_size).mean(axis=1) # Reshape the data into bins and calculate the mean of each bin

def calculate_binning_errors(data, bin_sizes, discard=0):
    """Estimate the standard error of the mean for different bin sizes."""
    data = np.asarray(data[discard:], dtype=float) # Discard the initial portion of the data before equilibration
    errors = []

    for bin_size in bin_sizes:
        num_bins = len(data) // bin_size # Calculate the number of complete bins that can be formed with the given bin size

        if num_bins < 2:
            errors.append(np.nan)
            continue

        bin_means = bin_data(data, bin_size) # Calculate the mean of each bin
        error = np.std(bin_means, ddof=1) / np.sqrt(len(bin_means)) # Calculate the standard error of the mean for the bin means
        errors.append(error)

    return np.asarray(errors) # Return the calculated errors

def analytical_mean_energy_per_spin(num_spins, temperature):
    """Return the analytical mean energy per spin."""
    return -(num_spins - 1) / num_spins * np.tanh(1 / temperature) # Equation derived in Derivations.ipynb


def analytical_heat_capacity_per_spin(num_spins, temperature):
    """Return the analytical heat capacity per spin."""
    return (num_spins - 1) / (num_spins * temperature**2) / np.cosh(1 / temperature)**2 # Equation derived in Derivations.ipynb


def analytical_mean_magnetization_squared(num_spins, temperature):
    """Return the analytical mean squared magnetization."""
    correlation_factor = np.tanh(1 / temperature)
    distances = np.arange(1, num_spins)
    return num_spins + 2 * np.sum((num_spins - distances) * correlation_factor**distances) # Equation derived in Derivations.ipynb


def analytical_susceptibility_per_spin(num_spins, temperature):
    """Return the analytical magnetic susceptibility per spin."""
    mean_magnetization_squared = analytical_mean_magnetization_squared(num_spins, temperature)
    return mean_magnetization_squared / (num_spins * temperature) # Equation derived in Derivations.ipynb