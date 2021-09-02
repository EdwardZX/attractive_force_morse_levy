import hoomd


hoomd.context.initialize("--mode=gpu")
from hoomd import mpcd

box = hoomd.data.boxdim(L=10.)
hoomd.init.read_snapshot(hoomd.data.make_snapshot(N=0,box=box))

system = mpcd.init.make_random(N=int(5*box.get_volume()), kT=1.0, seed=42)

mpcd.integrator(dt=0.1)
mpcd.stream.bulk(period=1)
# crashes
mpcd.collide.at(seed=42, period=1, kT=1.0)
# doesn't crash
#mpcd.collide.srd(seed=42, angle=130, period=1, kT=1.0)

hoomd.run(2000)
