from pathlib import Path
from typing import Callable

import ovito
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def cross_corr(x, y):
    """
    Calculates cross-correlation function of x and y using the 
    fast Fourier transform.
    :param x: array[float], data set 1
    :param y: array[float], data set 2
    :return: cf: array[float], cross-correlation function
    """
    N=len(x)
    F1 = np.fft.fft(x, n=2**(N*2 - 1).bit_length())
    F2 = np.fft.fft(y, n=2**(N*2 - 1).bit_length())
    PSD = F1 * F2.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   
    n=N*np.ones(N)-np.arange(0,N)
    cf = res/n
    return cf


def msd_fft_cross(r, k):
    """
    Calculates "MSD" (cross-correlations) using the fast Fourier transform.
    :param r: array[float], positions of atom type 1 over time
    :param k: array[float], positions of atom type 2 over time
    :return: msd: array[float], "MSD" over time
    """
    N=len(r)
    D=np.multiply(r,k).sum(axis=1)
    D=np.append(D,0)
    S2=sum([cross_corr(r[:, i], k[:,i]) for i in range(r.shape[1])])
    S3=sum([cross_corr(k[:, i], r[:,i]) for i in range(k.shape[1])])
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
        Q=Q-D[m-1]-D[N-m]
        S1[m]=Q/(N-m)
    msd = S1-S2-S3
    return msd


def total_displacement_modifier(frame: int, data: ovito.data.DataCollection) -> None:

    types = data.particles['Particle Type'][...]
    per_particle_displacements = data.particles['Displacement'][...]

    # get per-type data and reduce it
    for type_ in set(types):
        displacements_of_type = per_particle_displacements[types == type_]
        data.attributes[f'Total displacement {type_:.0f}'] = np.sum(displacements_of_type, axis=0)


def select_n_closest(particle_index: int, n: int) -> Callable:

    def wrapper(frame: int, data: ovito.data.DataCollection) -> None:

        selection = data.particles_.create_property('Selection', data=None)
        finder = ovito.data.NearestNeighborFinder(n, data)

        for neigh in finder.find(particle_index):
            selection[neigh.index] = 1

    return wrapper


def recenter_modifier(only_selected: bool = False) -> Callable:

    def wrapper(frame: int, data: ovito.data.DataCollection) -> None:
        
        if not only_selected:
            center_of_mass_displacement = np.mean(data.particles["Displacement"][...], axis=0)
        else:
            center_of_mass_displacement = np.mean(data.particles["Displacement"][...][data.particles.selection[...] == 1], axis=0)
        data.particles_.positions_ -= center_of_mass_displacement

    return wrapper


def vacancy():

    # timestep in ns
    timestep = 1.0e-6

    # load in modified run and recenter it
    vacancy_run = ovito.io.import_file(Path("lammps_data/vacancy_with_pd.dump"))
    vacancy_run.modifiers.append(ovito.modifiers.CalculateDisplacementsModifier())
    vacancy_run.modifiers.append(total_displacement_modifier)
    vacancy_run.modifiers.append(recenter_modifier(only_selected=False))

    # populate displacement and time arrays
    displacement = np.zeros((vacancy_run.source.num_frames, 3))
    time = np.zeros(vacancy_run.source.num_frames)
    for frame in range(vacancy_run.source.num_frames):
        
        data = vacancy_run.compute(frame)

        displacement[frame, :] = data.attributes["Total displacement 0"]
        time[frame] = timestep * data.attributes["Timestep"]

    # plot results
    msd = msd_fft_cross(displacement, displacement)
    plt.plot(time, msd)
    plt.xlabel("time")
    plt.ylabel("mean square displacement")
    plt.grid()
    plt.savefig("plots/vacancy.pdf")
    plt.close()


def interstitial():

    # timestep in ns
    timestep = 1.0e-6

    interstitial_run = ovito.io.import_file(Path("lammps_data/interstitial_with_pd.dump"))
    interstitial_run.modifiers.append(ovito.modifiers.CalculateDisplacementsModifier())
    interstitial_run.modifiers.append(total_displacement_modifier)

    # need to unselect the dumb-bell particles before recentering
    # select particles closest to defect and invert selection
    interstitial_run.modifiers.append(select_n_closest(particle_index=0, n=2))
    # invert selection and then recenter
    interstitial_run.modifiers.append(ovito.modifiers.InvertSelectionModifier())
    interstitial_run.modifiers.append(recenter_modifier(only_selected=True))

    displacement = np.zeros((interstitial_run.source.num_frames, 3))
    time = np.zeros(interstitial_run.source.num_frames)
    for frame in range(interstitial_run.source.num_frames):
        
        data = interstitial_run.compute(frame)
        displacement[frame, :] = data.attributes["Total displacement 0"]
        time[frame] = timestep * data.attributes["Timestep"]

    msd = msd_fft_cross(displacement, displacement)
    plt.plot(time, msd)
    plt.xlabel("time")
    plt.ylabel("mean square displacement")
    plt.grid()
    plt.savefig("plots/interstitial.pdf")
    plt.close()


if __name__ == "__main__":

    mpl.use("Agg")
    vacancy()
    interstitial()
