using Sunny, GLMakie, FFTW

latvecs = lattice_vectors(1, 1.1, 1.2, 90, 90, 90)
crystal = Crystal(latvecs, [[0,0,0]])
view_crystal(crystal; ndims=2)

################################################################################
# FM Interaction
################################################################################
L = 1
J = -1.0
sys = System(crystal, [1 => Moment(; s=1/2, g=2)], :dipole; dims=(L,1,1))
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)

swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
qs = q_space_path(crystal, [[-0.5, 0, 0], [0, 0, 0], [0.5, 0, 0]], 300)
is = intensities_bands(swt, qs)
plot_intensities(is)


################################################################################
# AFM Interaction
################################################################################
L = 2
J = 0.45
sys = System(crystal, [1 => Moment(; s=1/2, g=2)], :dipole; dims=(L,1,1))
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)

swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
qs = q_space_path(crystal, [[0, 1/2, 1/2], [1/4, 1/2, 1/2], [1/2, 1/2, 1/2], [3/4, 1/2, 1/2], [1, 1/2, 1/2]], 300)
energies = range(0, 1, 200)
is = intensities(swt, qs; energies, kernel=gaussian(; fwhm=0.05))
plot_intensities(is; colormap=:jet1)