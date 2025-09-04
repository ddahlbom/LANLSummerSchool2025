using DrWatson
@quickactivate "LANLSummerSchool2025"

using Sunny, GLMakie

crystal = Sunny.pyrochlore_crystal()
# view_crystal(crystal)

J1  = 38.05
J2  = 0.0815*J1
J3a = 0.1050*J1
J3b = 0.0085*J1 

dims = (10, 10, 4)
sys = System(crystal, [1 => Moment(s=1, g=2)], :dipole; dims)

set_exchange!(sys, J1,  Bond(1, 6, [0, 0, 0]))
set_exchange!(sys, J2, Bond(3, 5, [0, 0, 0]))
set_exchange!(sys, J3a, Bond(1, 3, [0, 0, 0]))
set_exchange!(sys, J3b, Bond(1, 3, [-1, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)
energy_per_site(sys)
plot_spins(sys)

################################################################################
# Static SF
################################################################################
kT = 0.01
damping = 0.1
dt = 0.001
integrator = Langevin(dt; damping, kT)

sssf = SampledCorrelationsStatic(sys; measure=ssf_perp(sys))
for _ in 1:5
    for _ in 1:100
        step!(sys, integrator)
    end
    minimize_energy!(sys)
    add_sample!(sssf, sys)
end

grid = q_space_grid(crystal, [1, 0, 0], range(-2.5, 2.5, 200), [0, 1, 0], (-2.5, 2.5))
res_sf = intensities_static(sssf, grid)
plot_intensities(res_sf)


################################################################################
# Dynamic SF
################################################################################
energies = range(0, 20, 201)
dssf = SampledCorrelations(sys; dt, energies, measure=ssf_trace(sys))

for _ in 1:100
    step!(sys, integrator)
end

@time add_sample!(dssf, sys)
path = q_space_path(crystal, [[-1/2, 0, 0], [0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]], 300)
res = intensities(dssf, path; energies=energies[1:end-1], kT)
plot_intensities(res)