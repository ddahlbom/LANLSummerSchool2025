using Sunny, GLMakie, FFTW

latvecs = lattice_vectors(1, 1, 1.2, 90, 90, 120)
crystal = Crystal(latvecs, [[0,0,0]])
view_crystal(crystal; ndims=2)

L = 21
sys = System(crystal, [1 => Moment(; s=1/2, g=-1)], :dipole; dims=(L, L, 1))
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))

################################################################################
# Animation
################################################################################
randomize_spins!(sys)
minimize_energy!(sys; maxiters=10_000)

fig = plot_spins(sys; colorfn=i->sys.dipoles[i][3], colorrange=(-0.9, 0.9), colormap=:roma, compass=false, show_cell=false)
lscene = fig.figure.content[1]
rotate_cam!(lscene.scene, (-π/4, 0, 0))
# translate_cam!(lscene.scene, Vec3f(-12, -0, 0))
zoom!(lscene.scene, 0.7)

for _ in 1:500
    step!(sys, integrator)
end

dt = 0.05
kT = 0.05
integrator = Langevin(dt; damping=0.1, kT)
mp_integrator = ImplicitMidpoint(dt)

# Perform the simulation and save the animation.
record(fig, "tl_animation.mp4", 1:270; framerate=30) do _
    notify(fig)
    for _ in 1:4
        step!(sys, mp_integrator)
    end
end


################################################################################
# Dynamic Structure Factor 
################################################################################
randomize_spins!(sys)
minimize_energy!(sys; maxiters=10_000)
for _ in 1:500
    step!(sys, integrator)
end

energies=range(0, 3.0, 200)
sc = SampledCorrelations(sys; measure=ssf_trace(sys), dt, energies)

for _ in 1:20
    for _ in 1:500
        step!(sys, integrator)
    end
    add_sample!(sc, sys)
end

path = q_space_path(crystal, [[-1/2, -1/2, 0], [0, 0, 0], [1/2, 1/2, 0], [1/2, 0, 0]], 400)
res = intensities(sc, path; energies, kT)
plot_intensities(res; axisopts=(; xticklabelrotation=π/4,))


################################################################################
# Static Structure Factor 
################################################################################
randomize_spins!(sys)
minimize_energy!(sys; maxiters=10_000)

sc = SampledCorrelationsStatic(sys; measure=ssf_trace(sys))
add_sample!(sc, sys)

grid = q_space_grid(crystal, [1, 0, 0], range(-1.0, 1.0, 300), [0, 1, 0], (-1.0, 1.0); orthogonalize=true)
res = intensities_static(sc, grid)
plot_intensities(res; axisopts=(; xticklabelrotation=π/4,), colorrange=(0, 100.0))