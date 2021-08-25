# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import hoomd.md
import numpy
from matplotlib import pyplot
import ex_render

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #print_hi('Hello hoomd PyCharm')
    '''
    data = numpy.genfromtxt(fname='log-output.log', skip_header=True)
    ex_render.display_movie(ex_render.render_sphere_frame, 'trajectory.gsd')
    '''


    hoomd.context.initialize('')
    hoomd.init.create_lattice(unitcell=hoomd.lattice.sc(a=2.0), n=5)
    nl = hoomd.md.nlist.cell()

    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
    lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)

    hoomd.md.integrate.mode_standard(dt=0.005)
    all = hoomd.group.all()



    hoomd.analyze.log(filename="log-output.log",
                      quantities=['potential_energy', 'temperature'],
                      period=100,
                      overwrite=True)

    hoomd.dump.gsd("trajectory.gsd", period=2e3, group=all, overwrite=True)
    hoomd.run(5e4)



# See PyCharm help at https://www.jetbrains.com/help/pycharm/

