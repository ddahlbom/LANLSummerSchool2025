using Sunny, GLMakie, FFTW

latvecs = lattice_vectors(1, 1.1, 1.2, 90, 90, 90)
crystal = Crystal(latvecs, [[0,0,0]])
view_crystal(crystal; ndims=2)

################################################################################
# FM Interaction
################################################################################
L = 10
J = -1.0
sys = System(crystal, [1 => Moment(; s=1/2, g=2)], :dipole; dims=(L,1,1))
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2, orthographic=true, show_cell=false)


################################################################################
# AFM Interaction
################################################################################
L = 10
J = 1.0
sys = System(crystal, [1 => Moment(; s=1, g=2)], :dipole; dims=(L,1,1))
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
# set_onsite_coupling!(sys, -0.1spin_matrices(1)[2]^2, 1)

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2, orthographic=true, show_cell=false)

################################################################################
# Single-ion interaction
################################################################################
L = 10
J = 1.0
sys = System(crystal, [1 => Moment(; s=1, g=2)], :dipole; dims=(L,1,1))
set_onsite_coupling!(sys, -0.1spin_matrices(1)[2]^2, 1)

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2, orthographic=true, show_cell=false)


################################################################################
# Competing NN and NNN
################################################################################
L = 40
sys = System(crystal, [1 => Moment(; s=1/2, g=2)], :dipole; dims=(L,1,1))
J1 = 1.0
J2 = 2J1
set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, J2, Bond(1, 1, [2, 0, 0]))

begin
    randomize_spins!(sys)
    minimize_energy!(sys)
    plot_spins(sys; ndims=2)
end

begin
    zs = [dipole[3] for dipole in sys.dipoles]
    zs_ft = fft(zs)
    sf = real.(zs_ft .* conj(zs_ft))[:] ./ prod(dims)


    qs = (0:L-1) / L
    xticks = ([0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0], ["0", "π/8", "π/4", "3π/8", "π/2", "5π/8", "3π/4", "7π/8", "2π"]) 
    fig = lines(qs, sf; axis=(xticks=xticks,))
    scatter!(qs, sf)
    fig
end


################################################################################
# Many competing interactions -- skyrmions
################################################################################
latvecs = lattice_vectors(1, 1, 1.1, 90, 90, 120)
crystal = Crystal(latvecs, [[0,0,0]])
view_crystal(crystal; ndims=2)
sys = System(crystal, [1 => Moment(; s=1, g=-1)], :dipole_large_S; dims=(12, 12, 1))
J1 = -1.0
J2 = 1/2
h = 0.15
D = 0.3/2
set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, J2, Bond(1, 1, [1, 2, 0]))
set_onsite_coupling!(sys, -D*spin_matrices(Inf)[3]^2, 1)
set_field!(sys, [0, 0, h])

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys)