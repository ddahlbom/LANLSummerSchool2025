using DrWatson
@quickactivate "LANLSummerSchool2025"

using Sunny, GLMakie

crystal = Sunny.pyrochlore_crystal()
view_crystal(crystal)

J1 = 3.2
J2 = 3.2
J3 = 0.019 
J4 = -0.07 
Jnn = [
     J2 J4 J4
    -J4 J1 J3
    -J4 J3 J1
]
# Jnn = [
#      J2 J3 J4
#      J3 J1 J4
#     -J4 -J4 J1
# ]
Jnnn = - 0.025


dims = (4, 4, 4)
sys = System(crystal, [1 => Moment(s=1, g=2)], :dipole; dims)

print_symmetry_table(crystal, 1.5)
set_exchange!(sys, Jnn,  Bond(6, 1, [0, 0, 0]))
set_exchange!(sys, Jnnn, Bond(3, 5, [0, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)
energy_per_site(sys)
plot_spins(sys)