import hoomd
import hoomd.md

import numpy as np
from scipy.stats import levy_stable,levy ## alpha ,beta means
import matplotlib.pyplot as plt

def get_ini_params(n,a):
    unitcell = hoomd.lattice.sq(a=a)
    snap = unitcell.get_snapshot()
    n = [n,n]
    if snap.box.dimensions == 3:
        snap.replicate(n[0], n[1], n[2])
    if snap.box.dimensions == 2:
        snap.replicate(n[0], n[1], 1)
    np.random.seed(123)
    particles_N = int(3 * snap.particles.N / 4)
    idx = np.random.choice(snap.particles.N, particles_N, replace=False)
    snap.particles.typeid[:] = 1
    snap.particles.typeid[idx] = 0
    snap.particles.types=['A','B']


    system = hoomd.init.read_snapshot(snapshot=snap)




    return particles_N, system

def my_init(init_t):
    all = hoomd.group.all()
    nl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2 ** (1 / 6), nlist=nl)
    lj.set_params(mode='shift')
    lj.pair_coeff.set(['A','B'], ['A','B'], epsilon=1.0, sigma=1.0)

    hoomd.md.integrate.mode_standard(dt=init_t)

    init_intgrate = hoomd.md.integrate.brownian(group=all, kT=0.85, seed=123)
    logger = hoomd.dump.gsd("init_setting.gsd", period=500, group=all, overwrite=True)
    hoomd.run(2500+1)
    return init_intgrate, logger

def my_update(mode,dt,num):
    nl = hoomd.md.nlist.cell()
    r_cut = 5
    group_A = hoomd.group.type('A')

    wca = hoomd.md.pair.lj(r_cut=2 ** (1 / 6), nlist=nl)
    wca.set_params(mode='shift')
    #wca.pair_coeff.set('A','A',epsilon=1.0, sigma=1.0)
    wca.pair_coeff.set(['A', 'B'], ['A', 'B'], epsilon=0.0, sigma=1.0)
    wca.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)

    # attractive morse chemistry
    morse = hoomd.md.pair.morse(r_cut=20, nlist=nl)
    morse.pair_coeff.set(['A', 'B'], ['A', 'B'], r_on=1)
    morse.pair_coeff.set('A', 'B', D0=2.5, alpha=0.25, r0=1.25)
    morse.pair_coeff.set('A', 'A', D0=0, alpha=1.0, r0=1.25)
    morse.pair_coeff.set('B', 'B', D0=0, alpha=1.0, r0=1.25)

    # physic attractive
    # nl.add_exclusions(hd.system.particles[0].tag, hd.system.particles[1].tag)
    # lj = hoomd.md.pair.lj(r_cut=r_cut, nlist=nl)
    # lj.pair_coeff.set(['A', 'B'], ['A', 'B'], r_on=1)
    # lj.pair_coeff.set('A', 'A', epsilon=0, sigma=2)
    # lj.pair_coeff.set('A', 'B', epsilon=0.31, sigma=2)
    # lj.pair_coeff.set('B', 'B', epsilon=0, sigma=0.1)
    #nl.reset_exclusions(exclusions=[])

    hoomd.md.integrate.mode_standard(dt=dt)
    hd_levy = hoomd.md.integrate.brownian(group=group_A, kT=0.25, seed=123)

    filename = ("trajectory_"
                + mode
                + ".gsd")
    logger = hoomd.dump.gsd(filename, period=2e3, group=hoomd.group.all(), overwrite=True)

    if mode == 'levy':
        print('levy mode:')
        tau = 300; force = None

        for i in range(0,int(1e5),tau):
            hd_levy.disable()
            np.random.seed(i)
            F_max = 1/0.0005
            F_norm = levy.rvs(size=num)
            F_norm[F_norm>F_max] = F_max
            theta = 2*np.pi * np.random.random(num)
            F_active = np.vstack([F_norm * np.cos(theta),F_norm*np.sin(theta)]).T
            zero_padding = np.zeros([num,1])
            F_active = list(map(tuple, 5e-3*np.hstack([F_active,zero_padding])))

            if force is not None:
                force.disable()
            force = hoomd.md.force.active(group=group_A,
                                  seed=123,
                                  f_lst=F_active,
                                  orientation_link=False)


            #hoomd.md.force._force

            #force.enable()
            #print('levy setting')

            # print(lv)
            # fig, ax = plt.subplots(1, 1)
            # ax.hist(lv, density=True, histtype='stepfilled', alpha=0.2)
            # plt.show()

            # hoomd.dump.gsd("trajectory_levy.gsd", period=tau-1, group=hoomd.group.all(), overwrite=False)
            # hoomd.md.integrate.mode_standard(dt=dt) # reset
            # hd_levy = hoomd.md.integrate.brownian(group=group_A, kT=0.25, seed=123)
            #force = hoomd.md.pair
            hd_levy.enable()
            hoomd.run(tau)

    else:
        #force = hoomd.md.pair
        hoomd.run(1e5)

    logger.disable()

    #hoomd.run(1e5)

if __name__ == '__main__':
    params = {'N':8,'a':15,'dt':0.001}
    hoomd.context.initialize("")

    particles_N,system = get_ini_params(params['N'],params['a'])
    init_intgrate, logger= my_init(params['dt'])
    init_intgrate.disable(); logger.disable()
    print('init finished')
    my_update(mode='all',dt=params['dt'],num=particles_N)
    print('finished')