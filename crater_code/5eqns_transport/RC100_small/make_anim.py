"""
100 m hemisphere
"""

import sys
if 'matplotlib' not in sys.modules:
    # Use an image backend to insure animation has size specified by figsize
    import matplotlib
    matplotlib.use('Agg')

from pylab import *
from matplotlib import animation
from scipy import io as sio
import os,sys,glob


root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'


def rel_import(fullpath):
    """
    Return a module imported from a full path name.
    To reload the module, call this again (rather than using importlib.reload).
    """

    import os, sys, importlib
    fname = os.path.split(fullpath)[1]
    modname = os.path.splitext(fname)[0]
    spec = importlib.util.spec_from_file_location(modname, fullpath)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[modname] = module
    print('loaded module from file: ',module.__file__)
    return module

C = rel_import(os.path.join(root_dir, 'compareALE3D.py'))

outdir_mfluid = os.path.join(root_dir, 'RC100_small/_output_5m_top2km')

datadir_ale = '/Users/rjl/D/Darrel_craters/surface_data_250804'
datadir_ale = '/Users/rjl/D/Darrel_craters/surface_data_250905'

fname_prefix_ale = 'hemi100m_4km'

tf_ale, find_frame_ale = C.load_times_ale(datadir_ale, fname_prefix_ale)
tf_mfluid, find_frame_mfluid = C.load_times_mfluid(outdir_mfluid)

# Create intial plot:

fig = figure(figsize=(12,7))

t = 0.


# mfluid:
frameno_mfluid = find_frame_mfluid(t)
rkm_mfluid, eta_mfluid, t_mfluid = C.load_surf_mfluid(outdir_mfluid,
                                                    frameno_mfluid)
mfluid_plot, = plot(rkm_mfluid, eta_mfluid, 'r', label='mfluid')


# ALE:
frameno_ale = find_frame_ale(t)
fname = '%s_surface_%s.mat' % (fname_prefix_ale, str(frameno_ale).zfill(4))
rkm, eta, t = C.load_surf_ale(datadir_ale, fname)
ale_plot, = plot(rkm, eta, 'b', label='ALE3D')


legend(loc='upper right', framealpha=1)
xlabel('radial distance (km)')
ylabel('surface elevation (m)')
xlim(0,1.5)
ylim(-150,200)
#ylim(-500,1000)
grid(True)
title_text = title('Hemispherical crater with radius 100 m on 4km ocean')


def update(t):
    """
    Update an existing plot by resetting the data.
    Use the frames from each simulation closest to the given time t.
    """

    # ALE:
    frameno_ale = find_frame_ale(t)
    fname = '%s_surface_%s.mat' % (fname_prefix_ale, str(frameno_ale).zfill(4))
    rkm, eta, t_ale = C.load_surf_ale(datadir_ale, fname)
    ale_plot.set_data(rkm, eta)

    # mfluid:
    frameno_mfluid = find_frame_mfluid(t)

    rkm, eta, t_mfluid = C.load_surf_mfluid(outdir_mfluid, frameno_mfluid)
    mfluid_plot.set_data(rkm, eta)
    if max(abs(t_ale-t),abs(t_mfluid-t)) > 0.1:
        # t not found exactly in one sim or the other:
        print('MISMATCH t = %.1f: frameno_ale = %i at t = %1.f, frameno_mfluid = %i at t = %1.f' \
            % (t,frameno_ale, t_ale, frameno_mfluid, t_mfluid))

    title_text.set_text('Hemispherical crater, radius 100 m, on 4km ocean,' \
            + '  t = %6.1f  (ALE3D at t = %6.1f)' % (t_mfluid, t_ale))


if __name__ == '__main__':

    print('Making anim...')
    times = tf_mfluid[:,1]
    #times = [0, 30, 60]
    anim = animation.FuncAnimation(fig, update, frames=times,
                                   interval=200, blit=False)

    # Output files:
    name = 'RC100_animation_newALE'

    fname_mp4 = name + '.mp4'

    fname_html = None
    #fname_html = name + '.html'

    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        animation_tools.make_html(anim, file_name=fname_html, title=name)
