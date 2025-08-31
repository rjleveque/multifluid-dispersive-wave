from pylab import *
import os,sys
from scipy.interpolate import interp1d

root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'


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

C = fullpath_import(os.path.join(root_dir, 'compareALE3D.py'))

AGT = os.environ['AGT']
LW = fullpath_import(os.path.join(AGT, 'linear_transforms', 'linear_waves.py'))


outdir_mfluid = 'RC1000_small/_output_40m_100km_t600'
h0 = 4000.
RC = 1000.
tAiry_start = 300.
tAiry_end = 4000.
rAiry_end = 100e3
savefile_name = 'eta_hu_bc_RC%i_h%s_tstart%i_rbc%ikm' \
    % (int(RC), int(h0), int(tAiry_start), int(rAiry_end/1000))
savefile_dir = './BC'
os.system('mkdir -p %s' % savefile_dir)

grav = 9.81

print('outdir_mfluid = ',outdir_mfluid)
tf_mfluid, find_frame_mfluid = C.load_times_mfluid(outdir_mfluid)

def load_mfluid(frameno):
    fname = outdir_mfluid + '/fort.c%s' % str(frameno).zfill(4)
    print('Loading ',fname)
    d = loadtxt(fname)
    r_mfluid = d[:,0]
    eta_mfluid = d[:,2] - h0
    u_mfluid = d[:,1]
    return r_mfluid, eta_mfluid, u_mfluid

def save_eta_u(t, r, eta, u, fname):
    d = vstack((r,eta,u)).T
    savetxt(fname, d, header='%.1f\n%i' % (t,len(eta)), comments='')
    print('Created ',fname)

xlower = 0
xupper = 200e3  # maximum distance to start of continental slope

# Initial eta and u
L = xupper - xlower
mx = 5000
dx = L/mx
r = linspace(dx/2, xupper-dx/2, mx)
k = linspace(1e-6,0.005,4000)

def make_Airy_data(tAiry_start=250, r=r, k=k):
    frameno_mfluid = find_frame_mfluid(tAiry_start)
    rkm_mfluid, eta_mfluid, u_mfluid, t_mfluid = \
                    C.load_mfluid(outdir_mfluid, frameno_mfluid)

    if t_mfluid != tAiry_start:
        print('Resetting tAiry_start from %s to %s' % (tAiry_start,t_mfluid))
        tAiry_start = t_mfluid
    r_mfluid = rkm_mfluid * 1000.

    if 0:
        # damp out near origin:
        eta_mfluid = where(r_mfluid<3e3,
                            eta_mfluid*exp(-((r_mfluid-3e3)/1e3)**2), eta_mfluid)
        u_mfluid = where(r_mfluid<3e3,
                            u_mfluid*exp(-((r_mfluid-3e3)/1e3)**2), u_mfluid)


    # resample eta_mfluid at r, possibly extending by 0 to larger distance:
    print('resampling eta_mfluid from r_mfluid (out to %s m with dr=%s)' \
            % (r_mfluid[-1], r_mfluid[-1]-r_mfluid[-2]) \
            + '\n                               to r (out to %s m with dr=%s)' \
            % (r[-1],dx))
    etafcn = interp1d(r_mfluid, eta_mfluid, fill_value=0., bounds_error=False)
    eta = etafcn(r)
    ufcn = interp1d(r_mfluid, u_mfluid, fill_value=0., bounds_error=False)
    u = ufcn(r)  # not needed?

    print('Computing Hankel transform...')
    etahat = LW.Htransform(r,eta,k)

    if 0:
        # NOT WORKING since etathat is complex
        d = vstack((k,etahat_mfluid)).T
        fname = 'etahat_mfluid_t%s.txt' % t
        savetxt(fname,d)
        print('Created ',fname)
    return eta, etahat, u

def load_etahat_mfluid(fname):
    # NOT WORKING since etathat is complex
    k,etahat_mfluid = loadtxt(fname, unpack=True)
    return k, etahat_mfluid

def evolve_Airy(etahat,tAiry_start,tAiry_end,k=k,r=r):
    # not used now, make time series at  rAiry_end instead for BC
    omega = lambda k: LW.omega_airy(k,h0)
    elapsed_time = tAiry_end - tAiry_start
    eta, u = LW.eta_u_radial(elapsed_time,r,k,etahat,omega,h0)
    return eta, u

# ========
# Create boundary data eta(rAiry_end,t) at fixed radius rAiry_end for series of times,
# and corresponding hu(rAiry_end,t)

def make_bc(rAiry_end, tAiry_end, tAiry_start, etahat_start, omega=None,
            save_txt=False):

    from scipy.signal import find_peaks
    dt = 1.
    t = arange(tAiry_start, tAiry_end+dt/2, dt)
    print('Computing time series...')
    elapsed_time = t - tAiry_start

    if omega is None:
        omega = lambda k: LW.omega_airy(k,h0)

    eta = LW.eta_tseries_radial(elapsed_time,rAiry_end,k,etahat_start,omega,h0)

    if 1:
        # for use with SWE on shelf:
        cphase = sqrt(grav*h0) # independent of k
        hu_swe = sqrt(grav*h0) * eta  # for SWE

        # for use with SGN on shelf:
        omega = lambda k: LW.omega_sgn(k,h0,alpha=1.153)
        jpeaks = find_peaks(eta)[0]
        jpeaks = jpeaks[4:]  # throw away first few points  - ADJUST?
        tpeaks = t[jpeaks]
        Tperiod = diff(t[jpeaks])
        tpeaks = hstack((0, tpeaks, tpeaks[-1]+1000))
        #Tperiod = hstack((240, Tperiod, Tperiod[-1], Tperiod[-1]))
        Tperiod = hstack((Tperiod[0], Tperiod, Tperiod[-1], Tperiod[-1]))
        Tfcn = interp1d(tpeaks, Tperiod) #, fill_value=0., bounds_error=False)
        Tt = Tfcn(t) # estimate of period at each time in t
        kk = linspace(0,0.01,1000)
        ww = omega(kk)
        kfcn = interp1d(ww,kk)  # k as a function of omega
        kt = kfcn(2*pi/Tt)  # estimate of k at each t
        cphase = omega(kt)/kt
        hu_sgn = cphase * eta
    else:
        raise NotImplementedError('have not implemented cgroup approach')
        # approximate group velocity of waves arriving at time t:
        #cgroup = (rAiry_end - 1) / t
        # need to invert for k and then compute cphase = omega(k)/k
        #hu = cphase * eta

    d = vstack((t,eta,hu_swe,hu_sgn)).T
    if save_txt:
        #fname = 'eta_hu_bc_%skm.txt' % int(rAiry_end/1e3)
        fname = os.path.join(savefile_dir, savefile_name + '.txt')
        savetxt(fname, d, header='%.0f\n%.1f\n%i' \
                                 % (rAiry_end,tAiry_start,len(eta)),
                comments='',fmt='%20.10e')
        print('Created ',fname)
    return t,eta,hu_swe,hu_sgn


def plot_bc(t,eta,hu_swe,hu_sgn,rAiry_end,tAiry_end,
            save_png=False):
    figure(figsize=(10,7))
    subplot(211)
    plot(t/60,eta,'b')
    title('time series of eta, hu at rAiry_end = %s km' \
            % int(rAiry_end/1e3), fontsize=15)
    grid(True)
    xlim(0,tAiry_end/60)
    #xlabel('time (minutes)', fontsize=13)
    ylabel('surface eta (m)', fontsize=13)
    xticks(fontsize=13)
    yticks(fontsize=13)

    subplot(212)
    plot(t/60,hu_swe,'b',label='hu for swe')
    plot(t/60,hu_sgn,'r',label='hu for sgn')
    #title('time series of hu(r2,t) at r2 = rAiry_end = %s km' \
    #        % int(rAiry_end/1e3), fontsize=15)
    grid(True)
    legend(loc='upper left', framealpha=1, fontsize=11)
    xlim(0,tAiry_end/60)
    xlabel('time (minutes)', fontsize=13)
    ylabel('momentum hu (m**2/s)', fontsize=13)
    xticks(fontsize=13)
    yticks(fontsize=13)

    tight_layout()

    if save_png:
        #fname = 'eta_hu_bc_%skm.png' % int(rAiry_end/1e3)
        fname = os.path.join(savefile_dir, savefile_name + '.png')
        savefig(fname, bbox_inches='tight')
        print('Created ',fname)

def test_bc(tAiry_start=250, etahat_start=None):

    if etahat_start is None:
        eta_start, etahat_start, u_start = make_Airy_data(tAiry_start=tAiry_start,
                                                          r=r, k=k)
    #rAiry_end = 200.e3
    #tAiry_end = 3.*3600.
    rAiry_end = 100.e3
    tAiry_end = 4000.

    t,eta,hu_swe,hu_sgn = make_bc(rAiry_end, tAiry_end, tAiry_start, etahat_start,
                       save_txt=True)
    plot_bc(t, eta,hu_swe,hu_sgn, rAiry_end, tAiry_end, save_png=True)

    return t,eta,hu_swe,hu_sgn

def compare_airy_sgn(etahat_start=None):

    tAiry_start = 250
    if etahat_start is None:
        eta_start, etahat_start, u_start = make_Airy_data(tAiry_start=tAiry_start,
                                                          r=r, k=k)
    rAiry_end = 200.e3
    tAiry_end = 3*3600.
    print('tAiry_end = %.1f sec' % tAiry_end)

    omega = lambda k: LW.omega_sgn(k,h0,1.153)
    t,eta,hu = make_bc(rAiry_end, tAiry_end, tAiry_start,
                       etahat_start,omega=omega, save_txt=False)
    plot_bc(t, eta, rAiry_end, tAiry_end, plot_eta=False)
    plot(t/60,eta,'r',label='SGN')

    t,eta,hu = make_bc(rAiry_end, tAiry_end, tAiry_start, etahat_start,
                       save_txt=False)
    #plot_bc(t, eta, rAiry_end, tAiry_end, plot_eta=False)
    plot(t/60,eta,'b',label='Airy')
    legend(loc='upper right', framealpha=1, fontsize=12)

    return t,eta,hu

if __name__ == '__main__':

    eta_start, etahat_start, u_start = make_Airy_data(tAiry_start=tAiry_start,
                                                      r=r, k=k)

    t,eta,hu_swe,hu_sgn = make_bc(rAiry_end, tAiry_end, tAiry_start,
                                  etahat_start, save_txt=True)

    plot_bc(t, eta,hu_swe,hu_sgn, rAiry_end, tAiry_end, save_png=True)
