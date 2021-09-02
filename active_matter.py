import hoomd
import hoomd.md


hoomd.context.initialize("--mode=gpu")


system = hoomd.init.create_lattice(unitcell=hoomd.lattice.sq(a=5),
                                   n=15)
all = hoomd.group.all()
N = len(all)
nl = hoomd.md.nlist.cell();
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl);
lj.set_params(mode='shift');

lj.pair_coeff.set('A', 'A', epsilon=1, sigma=1);

import numpy
activity = [ ( ((numpy.random.rand(1)[0] - 0.5) * 2),
               ((numpy.random.rand(1)[0] - 0.5) * 2),
               0)
             for i in range(N)];



hoomd.md.force.active(group=all,
                      seed=123,
                      f_lst=activity,
                      rotation_diff=0.005,
                      orientation_link=False);

hoomd.md.integrate.mode_standard(dt=1e-3);

hoomd.md.integrate.brownian(group=all, kT=0.2, seed=20);
hoomd.dump.gsd("trajectory.gsd", period=500, group=all, overwrite=True);
hoomd.run(1);