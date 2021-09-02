import hoomd
import hoomd.md


hoomd.context.initialize("--mode=cpu")


snapshot = hoomd.data.make_snapshot(N=10,
                                    #box=hoomd.data.boxdim(Lx=10, Ly=0.5, Lz=0.5),
                                    box=hoomd.data.boxdim(Lx=10, Ly=0.5, Lz=0.5),
                                    particle_types=['A', 'B'],
                                    bond_types=['polymer'])


snapshot.particles.position[:] = [[-4.5, 0, 0], [-3.5, 0, 0],
                                  [-2.5, 0, 0], [-1.5, 0, 0],
                                  [-0.5, 0, 0], [0.5, 0, 0],
                                  [1.5, 0, 0], [2.5, 0, 0],
                                  [3.5, 0, 0], [4.5, 0, 0]]


snapshot.particles.typeid[0:7]=0
snapshot.particles.typeid[7:10]=1

snapshot.bonds.resize(9)
snapshot.bonds.group[:] = [[0,1], [1, 2], [2,3],
                           [3,4], [4,5], [5,6],
                           [6,7], [7,8], [8,9]]


snapshot.replicate(1,20,20)
#snapshot.box = hoomd.data.boxdim(Lx=20, Ly=20, Lz=20)
hoomd.init.read_snapshot(snapshot)

nl = hoomd.md.nlist.cell()
dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, kT=0.8, seed=1)

dpd.pair_coeff.set('A', 'A', A=25.0, gamma = 1.0)
dpd.pair_coeff.set('A', 'B', A=100.0, gamma = 1.0)
dpd.pair_coeff.set('B', 'B', A=25.0, gamma = 1.0)

nl.reset_exclusions(exclusions = [])
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=100.0, r0=0)


hoomd.md.integrate.mode_standard(dt=1e-2)
all = hoomd.group.all()
integrator = hoomd.md.integrate.nve(group=all)#hoomd.group.type('B'))

integrator.randomize_velocities(kT=0.8, seed=42)

hoomd.analyze.log(filename="log-output.log",
                  quantities=['potential_energy', 'temperature'],
                  period=500,
                  overwrite=True)


hoomd.dump.gsd("trajectory.gsd", period=10e2, group=all, overwrite=True)

hoomd.run(5e3)

#import ex_render
#ex_render.display_movie(ex_render.render_sphere_frame, 'trajectory.gsd')