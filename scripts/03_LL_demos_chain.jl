using Sunny, GLMakie, FFTW

latvecs = lattice_vectors(1, 1.1, 1.2, 90, 90, 90)
crystal = Crystal(latvecs, [[0,0,0]])
view_crystal(crystal; ndims=2)

################################################################################
# De/coupled Interaction
################################################################################
L = 20
J = -1.0
sys = System(crystal, [1 => Moment(; s=1/2, g=-1)], :dipole; dims=(L,1,1))
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
set_field!(sys, [0, 0, 1.0])

dt = 0.01
kT = 0.2
integrator = Langevin(dt; damping=0.1, kT)
integrator_md = ImplicitMidpoint(dt)

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys)

for _ in 1:100
    step!(sys, integrator)
end
plot_spins(sys)


fig = plot_spins(sys; colorfn=i->sys.dipoles[i][3], colorrange=(-1, 1), ndims=2)
for _ in 1:500
    for _ in 1:5
        step!(sys, integrator_md)
    end
    notify(fig)
    sleep(1/60)
end

energies = range(0, 3.0, 200)
dssf = SampledCorrelations(sys; measure=ssf_trace(sys), dt, energies)

for _ in 1:20
    for _ in 1:100
        step!(sys, integrator)
    end
    add_sample!(dssf, sys)
end

path = q_space_path(crystal, [[-1/2, 0, 0], [0, 0, 0], [1/2, 0, 0]], 300)
res = intensities(dssf, path; energies, kT)
plot_intensities(res)






################################################################################
# Spin Wave Calc
################################################################################


swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
qs = q_space_path(crystal, [[-0.5, 0, 0], [0, 0, 0], [0.5, 0, 0]], 300)
is = intensities_bands(swt, qs)
plot_intensities(is)
