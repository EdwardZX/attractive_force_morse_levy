# 2 dimension of system of DNA-base and Au-particles;
# regular or random postion
# force:morse and lj(vander val)//constraint planar
# method: brownian update
import poisson_disc
import hoomd
import hoomd.md
import fresnel
import numpy as np
import ex_render
import matplotlib.pyplot as plt



if __name__ == '__main__':

    hoomd.context.initialize("")
    base_params = dict()
    np.random.seed(1) # fixed num

    ex_render.display_movie(ex_render.render_disk_frame, 'trajectory.gsd')
    L = 50
    r = 0.5
    space = 1.4
    base_params['rd'] = 0.1
    #N_base = m


    #
    #x = np.linspace(-L / 2, L / 2, m, endpoint=False)
    #position = list(itertools.product(x, repeat=3))
    # the position of base
    #pos = L*np.random.rand(base_params['N'],10)

    pos = 0.95*L * poisson_disc.Bridson_sampling(dims=np.array([1.0,1.0]),
                                            radius=5*7*2*base_params['rd']/L) - np.array([L/2,L/2])
    plt.scatter(pos[:,0],pos[:,1],s=20)
    plt.show()
    base_params['N'] = len(pos)
    pos_z = np.zeros([base_params['N'],1])
    base_params['pos'] = tuple(map(tuple, np.hstack([pos,pos_z])))
    #base_params['pos'] = tuple(map(tuple, pos))


    snapshot = hoomd.data.make_snapshot(N = base_params['N'],
                                        box=hoomd.data.boxdim(Lx=L,Ly=L,dimensions=2),
                                        particle_types=['A','B'],
                                        )

    snapshot.particles.position[:] = base_params['pos']

    idx = np.random.choice(base_params['N'], int(3*base_params['N']/4), replace=False)

    snapshot.particles.typeid[:] = 1
    snapshot.particles.typeid[idx] = 0
    #snapshot.particles.typeid[int(base_params['N']/4):] = 1 # which sorted by position
    snapshot.particles.diameter[:] = 2*base_params['rd'] * snapshot.particles.diameter


    hd = hoomd.init.read_snapshot(snapshot)


    #define the force
    nl = hoomd.md.nlist.cell()
    #nl.add_exclusions(hd.system.particles[0].tag, hd.system.particles[1].tag)
    lj = hoomd.md.pair.lj(r_cut=1.0,nlist=nl)
    lj.pair_coeff.set(['A','B'],['A','B'],r_on = r)
    lj.pair_coeff.set('A','A',epsilon=0.003974, sigma=0.1)
    lj.pair_coeff.set('A','B',epsilon=0.3974, sigma=0.1)
    lj.pair_coeff.set('B', 'B', epsilon=0, sigma=0.1)
    #hoomd.md.constrain.

    morse = hoomd.md.pair.morse(r_cut=1.0,nlist=nl)
    morse.pair_coeff.set(['A','B'],['A','B'],r_on = r)
    morse.pair_coeff.set('A', 'B',D0=1.0, alpha=1.0, r0=1.25*r)
    morse.pair_coeff.set('A', 'A', D0=0.001, alpha=1.0, r0=1.25*r,r_cut = 0)
    morse.pair_coeff.set('B', 'B', D0=0.001, alpha=1.0, r0=1.25*r,r_cut = 0)

    nl.reset_exclusions(exclusions=[])

    ## set the constrain
    groupA = hoomd.group.type('A')

    hoomd.md.integrate.mode_standard(dt=1e-3)

    integrator = hoomd.md.integrate.brownian(group=groupA,kT=0.1,dscale=2,seed=1)


    hoomd.dump.gsd("trajectory.gsd", period=5e2, group=hoomd.group.all(), overwrite=True)

    hoomd.run(1e4)

    #rigid = hoomd.md.constrain.rigid()

    print('finished')
