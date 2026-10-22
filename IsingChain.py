import numpy as np

class IsingChain:
    """Class representing a one-dimensional Ising model (chain) with open boundary conditions in the absence of an external magnetic field. J=1, k_B=1"""
    def __init__(self, num_spins, temperature, initial_chain=None):
        if num_spins < 1:
            raise ValueError("num_spins must be positive.")
        if temperature <= 0:
            raise ValueError("temperature must be positive.")

        self.num_spins = num_spins
        self.temperature = temperature

        if initial_chain is None:
            self.chain = np.random.choice([-1, 1], size=num_spins) # Choose random spins if no initial chain is provided
        else:
            self.chain = self.validate_chain(initial_chain) # Set chain to initial chain if provided, but validate input first

        self.energy = -np.sum(self.chain[:-1] * self.chain[1:]) # Initial energy
        self.magnetization = np.sum(self.chain) # Initial magnetization

        self.clear_history() # Set up history

    def validate_chain(self, chain):
        """Validate and return a copy of a spin chain."""
        chain = np.asarray(chain, dtype=int).copy()

        if chain.ndim != 1 or len(chain) != self.num_spins:
            raise ValueError("Chain must have length num_spins.")

        if not np.all(np.isin(chain, [-1, 1])):
            raise ValueError("In an Ising chain, every spin must be either -1 or +1.")

        return chain

    def clear_history(self):
        """Clear history and record the current state as the new initial state."""
        self.energy_history = [self.energy]
        self.magnetization_history = [self.magnetization]

    def change_chain(self, new_chain):
        """Replace the spin chain and reset its history."""
        self.chain = self.validate_chain(new_chain)

        self.energy = -np.sum(self.chain[:-1] * self.chain[1:])
        self.magnetization = np.sum(self.chain)

        self.clear_history()

    def energy_difference(self, index):
        """Calculate the energy change from flipping a spin."""
        left = self.chain[index - 1] if index > 0 else 0 # Open boundary conditions
        right = self.chain[index + 1] if index < self.num_spins - 1 else 0 # Open boundary conditions

        return 2 * self.chain[index] * (left + right) # This formula is derived in Derivations.ipynb

    def metropolis_step(self):
        """Perform a single Metropolis step."""
        index = np.random.randint(self.num_spins) # Find random spin to flip
        delta_energy = self.energy_difference(index) # Calculate energy change

        if delta_energy <= 0 or np.random.rand() < np.exp(-delta_energy / self.temperature): # Always accept ernergy-lowering flips. Accept energy-increasing flips with the Boltzmann probability exp(-delta_energy / temperature)
            self.magnetization -= 2 * self.chain[index] # Change magnetization
            self.energy += delta_energy # Change energy
            self.chain[index] *= -1 # Flip the spin

        # Append step to history
        self.energy_history.append(self.energy)
        self.magnetization_history.append(self.magnetization)

    def simulate(self, num_steps):
        """Run num_steps Metropolis steps and append to history."""
        if num_steps < 0:
            raise ValueError("num_steps must be non-negative.")

        for _ in range(num_steps):
            self.metropolis_step()

        return self.get_history()

    def get_history(self):
        """Return energy and magnetization histories as numpy arrays."""
        return np.asarray(self.energy_history), np.asarray(self.magnetization_history)

    def calculate_observables(self, discard=0):
        """Calculate observables from the recorded simulation history. Return a dictionary with 
        mean_energy,
        mean_energy_per_spin,
        mean_magnetization,
        mean_magnetization_per_spin,
        mean_magnetization_squared,
        mean_magnetization_squared_per_spin_squared,
        energy_variance,
        magnetization_variance,
        heat_capacity_per_spin,
        susceptibility_per_spin."""
        if discard < 0 or discard >= len(self.energy_history):
            raise ValueError("discard must be non-negative and smaller than the history length.")

        energy = np.asarray(self.energy_history[discard:], dtype=float)
        magnetization = np.asarray(self.magnetization_history[discard:], dtype=float)

        mean_energy = np.mean(energy)
        mean_magnetization = np.mean(magnetization)

        energy_variance = np.var(energy)
        magnetization_variance = np.var(magnetization)

        return {
            "mean_energy": mean_energy,
            "mean_energy_per_spin": mean_energy / self.num_spins,
            "mean_magnetization": mean_magnetization,
            "mean_magnetization_per_spin": mean_magnetization / self.num_spins,
            "mean_magnetization_squared": np.mean(magnetization**2),
            "mean_magnetization_squared_per_spin_squared": np.mean(magnetization**2) / self.num_spins**2,
            "energy_variance": energy_variance,
            "magnetization_variance": magnetization_variance,
            "heat_capacity_per_spin": energy_variance / (self.num_spins * self.temperature**2), # This formula is derived in Derivations.ipynb
            "susceptibility_per_spin": magnetization_variance / (self.num_spins * self.temperature), # This formula is derived in Derivations.ipynb
        }

