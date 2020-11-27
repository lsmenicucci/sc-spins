# import packages
import numpy as np


class TemperingHandler:
    def __init__(self, spintronics, betas):
        n = len(spintronics.rx)

        self.configurations = np.ones([betas, 3, n])
        self.buffer = np.zeros([3, n])
        self.spintronics = spintronics

    def load_configuration(self, idx):
        self.spintronics.rx = self.configurations[idx, 0, :]
        self.spintronics.ry = self.configurations[idx, 1, :]
        self.spintronics.rz = self.configurations[idx, 2, :]

    def unload_configuration(self, idx):
        self.configurations[idx, 0, :] = self.spintronics.rx
        self.configurations[idx, 1, :] = self.spintronics.ry
        self.configurations[idx, 2, :] = self.spintronics.rz

    def energy(self, idx):
        self.load_configuration(idx)
        return self.spintronics.total_energy()

    def swap(self, one, other):
        self.buffer = self.configurations[one, :]
        self.configurations[one, :] = self.configurations[other, :]
        self.configurations[other, :] = self.buffer

    def tempering_swap(self, one, other):
        one_energy = self.energy(one)
        other_energy = self.other(other)

        if(other_energy < one_energy):
            self.swap(one, other)
        else:
            r = np.random.random()
            d_beta = self.betas[one] - self.betas[other]

            if(r < np.exp((one_energy - other_enegy)*d_beta)):
                self.swap(one, other)

    def single_metropolis(self, idx, nsteps):
        self.load_configuration(idx)
        self.spintronics.metropolis(self.betas[idx], nsteps)
        self.unload_configuration(idx)

    def tempering_metropolis(self, nsteps):
        measures = []
        half_nsteps = int(nsteps/2.0)

        for idx in range(len(self.betas)):
            self.single_metropolis(self, idx, half_nsteps)

            if(idx > 0):
                self.single_metropolis(self, idx, half_nsteps)

                self.tempering_swap(idx - 1, idx)

                # Finish measurements for the previus idx
                self.single_metropolis(self, idx - 1, half_nsteps)

                # Gather measures
                measures_zip = zip(["mean_energy", "mean_mag_x", "mean_mag_y", "mean_mag_z"],
                                   spintronics.metropolis_measures)
                metropolis_measures = {n: v for (n, v) in measures_zip}
                metropolis_measures['beta'] = self.betas[idx - 1]
                measures.append(metropolis_measures)

            # Finish measurements for last
            if(idx == len(self.betas) - 1):
                self.single_metropolis(self, idx, half_nsteps)

                # Gather measures
                measures_zip = zip(["mean_energy", "mean_mag_x", "mean_mag_y", "mean_mag_z"],
                                   spintronics.metropolis_measures)
                metropolis_measures = {n: v for (n, v) in measures_zip}
                metropolis_measures['beta'] = self.betas[idx]
                measures.append(metropolis_measures)

        return measures
