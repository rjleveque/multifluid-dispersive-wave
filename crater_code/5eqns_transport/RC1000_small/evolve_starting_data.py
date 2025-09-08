import sys
if 'matplotlib' not in sys.modules:
    # Use an image backend to insure animation has size specified by figsize
    import matplotlib
    matplotlib.use('Agg')

from pylab import *
from matplotlib import animation

import os,sys
#import linear_waves as LW
from scipy.interpolate import interp1d

def fullpath_import(fullpath):
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

AGT = os.environ['AGT']
LW = fullpath_import(os.path.join(AGT, 'linear_transforms','linear_waves.py'))
C = fullpath_import(os.path.abspath('../compareALE3D.py'))

xlower = 0; xupper = 100e3; h0 = 4000.
# Initial eta and u
L = xupper - xlower
mx_H1 = 12000
dx = L/mx_H1
r = linspace(dx/2, xupper-dx/2, mx_H1)
k = linspace(1e-6,0.005,4000)

if 0:
    # wave packet:
    width = 10e3; r0 = 40e3; wavelength = 5e3; ampl = 100.
    eta0 = ampl*exp(-((r-r0)/width)**2) * cos(r*2*pi/wavelength)

#starting_data = loadtxt('starting.data',skiprows=2)

outdir = '_output_40m_100km_t600'
tf_mfluid, find_frame_mfluid = C.load_times_mfluid(outdir)

# initial time for Airy:
t0airy = 300
j = find_frame_mfluid(time=t0airy)
rkm, eta0, t0 = C.load_surf_mfluid(outdir, j)

r = rkm  * 1000.

etafcn = interp1d(r, eta0, fill_value=0., bounds_error=False)
#plot(r/1e3,eta0,'k--',label='orig eta at t0=%.0f' % t0)

if 0:
    # damp out near origin and for large r:
    eta0 = etafcn(r)
    #eta0 = where(r>35e3, eta0*exp(-(r-35e3)/5e3), eta0)
    eta0 = where(r<4e3, eta0*exp((r-4e3)/1e3), eta0)

print('Computing eta0hat transform...')
eta0hat = LW.Htransform(r,eta0,k)
omega = lambda k: LW.omega_airy(k,h0)
reval = r

if 0:
    figure(3, figsize=(9,5))
    clf()
    plot(k, eta0hat, 'k')
    xlabel('wavenumber k')
    ylabel('eta0hat')
    grid(True)
    title('Hankel transform of eta0')


# Create intial plot:

fig = figure(figsize=(12,7))

t = 0.

# mfluid:
frameno_mfluid = find_frame_mfluid(t)
rkm_mfluid, eta_mfluid, t_mfluid = C.load_surf_mfluid(outdir,
                                                    frameno_mfluid)
eta_airy = nan*eta_mfluid
airy_plot, = plot(rkm_mfluid, eta_airy, 'b',
         label='Airy, starting from mfluid at t=%.0f' % t0airy)
mfluid_plot, = plot(rkm_mfluid, eta_mfluid, 'r',
                    label='mfluid')
grid(True)
xlabel('distance from crater (km)')
ylabel('surface elevation (m)')
xlim(0,75.)
#ylim(-1200,1200)
ylim(-100,100)
grid(True)
title_text = title('Hemispherical crater, radius 1000 m, on 4km ocean,' \
            + '  t = %6.1f' % t_mfluid)

if 0:
    print('Evaluating eta(r,t) with Airy...')
    t = 600.
    omega = lambda k: LW.omega_airy(k,h0)

    eta,u = LW.eta_u_radial(t-t0,reval,k,eta0hat,omega,h0)
    plot(reval/1e3, eta, 'b', label='Airy evolution to t=%.0f' % t)

    j = find_frame(time=t)
    rkm, eta_mfluid, t_mfluid = C.load_surf_mfluid(outdir, j)
    airy_plot, = plot(rkm, eta_mfluid, 'b',
         label='Airy at t=%.0f starting from mfluid at t=%.0f' \
                % (t_mfluid,t0airy))
    mfluid_plot, = plot(rkm, eta_mfluid, 'r', label='mfluid at t=%.0f' % t_mfluid)

legend(loc='upper right', framealpha=1)


def update(t):
    """
    Update an existing plot by resetting the data.
    Use the frames from each simulation closest to the given time t.
    """

    # mfluid:
    frameno_mfluid = find_frame_mfluid(t)

    rkm, eta, t_mfluid = C.load_surf_mfluid(outdir, frameno_mfluid)
    mfluid_plot.set_data(rkm, eta)

    title_text.set_text('Hemispherical crater, radius 1000 m, on 4km ocean,' \
            + '  t = %6.1f' % t_mfluid)


    if t > t0airy:
        print('Evaluating eta(r,t=%.0f) with Airy...' % t)
        omega = lambda k: LW.omega_airy(k,h0)

        eta_airy,u_airy = LW.eta_u_radial(t-t0,reval,k,eta0hat,omega,h0)
        airy_plot.set_data(reval/1e3, eta_airy)


if __name__ == '__main__':

    print('Making anim...')
    times = tf_mfluid[:,1]
    #times = [0, 30, 60]
    anim = animation.FuncAnimation(fig, update, frames=times,
                                   interval=200, blit=False)

    # Output files:
    name = 'RC1000_AirySwitch_t%s' % str(t0airy).zfill(3)

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
