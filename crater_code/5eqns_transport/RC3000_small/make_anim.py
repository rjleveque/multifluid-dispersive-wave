

import sys
if 'matplotlib' not in sys.modules:
    # Use an image backend to insure animation has size specified by figsize
    import matplotlib
    matplotlib.use('Agg')

from pylab import *
from matplotlib import animation
from scipy import io as sio
import os,sys,glob


def load_surf_ale(datadir_ale, fname):
    """
    Load surface from Darrel's files, ALE3D results
    """
    mat_data = sio.loadmat(os.path.join(datadir_ale, fname))
    r = squeeze(array(mat_data['xs']))
    rkm = r / 1e3
    eta = squeeze(array(mat_data['ys']))
    try:
        time = float(mat_data['time'][0][0])
    except:
        print('time not found in .mat file')
        time = nan
    return rkm, eta, time


def load_times_ale(datadir_ale, fname_prefix_ale):
    files = glob.glob('%s/%s*' % (datadir_ale, fname_prefix_ale))
    files.sort()
    times = []

    for f in files:
        mat_data = sio.loadmat(os.path.join(datadir_ale, f))
        t = float(mat_data['time'][0][0])
        frameno = int(f[-8:-4])
        #print('ALE Frame %s: t = %.2f' % (frameno,t))
        times.append((frameno,float(round(t,1))))

    if 0:
        print('(frameno,time) for ALE results:')
        print(times)

    tf = array(times)

    def find_frame(time):
        """
        find frameno for best match to time
        """
        if time < tf[0,1]:
            k = 0
        else:
            k = where(tf[:,1] < time+1e-6)[0].max()
            #print('+++ k = %i' % k, '   tf[k:k+2, :] = ',tf[k:k+2, :])
            if (k < tf.shape[0]-1):
                if (tf[k+1,1] - time) < (time - tf[k,1]):
                    k = k+1
        return int(tf[k,0])

    return tf, find_frame

def load_surf_mfluid(outdir, j):
    h0 = 4000.
    fortc = '%s/fort.c%s' % (outdir, str(j).zfill(4))
    d = loadtxt(fortc)
    rkm = d[:,0] / 1e3
    eta = d[:,2] - h0
    fortt = fortc.replace('fort.c','fort.t')
    time = float(open(fortt).readline().split()[0])
    #print('frame %i: loaded fort.c at time %.1f sec' % (j,time))
    return rkm, eta, time

def load_times_mfluid(outdir):
    fortt_files = glob.glob('%s/fort.t*' % outdir)
    fortt_files.sort()
    #print(fortt_files)
    times = []
    for f in fortt_files:
        frameno = int(f[-4:])
        t = float(open(f).readline().split()[0])
        times.append((frameno,float(round(t,1))))
    tf = array(times)

    def find_frame(time):
        """
        find frameno for best match to time
        """
        if time < tf[0,1]:
            k = 0
        else:
            k = where(tf[:,1] < time+1e-6)[0].max()
            #print('+++ k = %i' % k, '   tf[k:k+2, :] = ',tf[k:k+2, :])
            if (k < tf.shape[0]-1):
                if (tf[k+1,1] - time) < (time - tf[k,1]):
                    k = k+1
        return int(tf[k,0])

    return tf, find_frame

def make_fig(title,xlim=(0,30.),ylim=(-1000,1000)):
    figure(figsize=(12,7))
    xlabel('radial distance (km)')
    ylabel('surface elevation (m)')
    xlim(0,30.)
    ylim(-1000,3000)
    grid(True)
    title(title_)

def compare_RC3000():

    outdir_mfluid = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/' \
                        + 'RC3000_small/_output_80m_32km'

    datadir_ale = '/Users/rjl/D/Darrel_craters/surface_data_250804'
    fname_prefix_ale = 'hemi3000m_4km'

    tf_ale, find_frame_ale = load_times_ale(datadir_ale, fname_prefix_ale)
    tf_mfluid, find_frame_mfluid = load_times_mfluid(outdir_mfluid)

    for k,t in enumerate([0,30,60]):
        figure(k, figsize=(12,7))
        clf()

        # ALE:
        frameno_ale = find_frame_ale(t)
        fname = '%s_surface_%s.mat' % (fname_prefix_ale, str(frameno_ale).zfill(4))
        rkm, eta, t = load_surf_ale(datadir_ale, fname)
        plot(rkm, eta, 'b', label='ALE frame %i at t = %.1f' % (frameno_ale,t))

        # mfluid:
        frameno_mfluid = find_frame_mfluid(t)
        rkm_mfluid, eta_mfluid, t_mfluid = load_surf_mfluid(outdir_mfluid,
                                                            frameno_mfluid)
        plot(rkm_mfluid, eta_mfluid, 'r', label='mfluid frame %i at t = %.1f' \
                % (frameno_mfluid,t_mfluid))

        legend(loc='upper right', framealpha=1)
        xlabel('radial distance (km)')
        ylabel('surface elevation (m)')
        xlim(0,30.)
        #ylim(-3000,3000)
        grid(True)
        title('Hemispherical crater with radius 300 m on 4km ocean')


if __name__ == '__main__':

    outdir_mfluid = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/' \
                        + 'RC3000_small/_output_80m_50km'

    datadir_ale = '/Users/rjl/D/Darrel_craters/surface_data_250804'
    fname_prefix_ale = 'hemi3000m_4km'

    tf_ale, find_frame_ale = load_times_ale(datadir_ale, fname_prefix_ale)
    tf_mfluid, find_frame_mfluid = load_times_mfluid(outdir_mfluid)

    fig = figure(figsize=(12,7))

    t = 0.

    # ALE:
    frameno_ale = find_frame_ale(t)
    fname = '%s_surface_%s.mat' % (fname_prefix_ale, str(frameno_ale).zfill(4))
    rkm, eta, t = load_surf_ale(datadir_ale, fname)
    ale_plot, = plot(rkm, eta, 'b', label='ALE3D')

    # mfluid:
    frameno_mfluid = find_frame_mfluid(t)
    rkm_mfluid, eta_mfluid, t_mfluid = load_surf_mfluid(outdir_mfluid,
                                                        frameno_mfluid)
    mfluid_plot, = plot(rkm_mfluid, eta_mfluid, 'r', label='mfluid')

    legend(loc='upper right', framealpha=1)
    xlabel('radial distance (km)')
    ylabel('surface elevation (m)')
    xlim(0,50.)
    #ylim(-3500,6000)
    ylim(-500,1000)
    grid(True)
    title_text = title('Hemispherical crater with radius 3000 m on 4km ocean')

    def update(t):
        """
        Update an existing plot
        """

        # ALE:
        frameno_ale = find_frame_ale(t)
        fname = '%s_surface_%s.mat' % (fname_prefix_ale, str(frameno_ale).zfill(4))
        rkm, eta, t_ale = load_surf_ale(datadir_ale, fname)
        ale_plot.set_data(rkm, eta)

        # mfluid:
        frameno_mfluid = find_frame_mfluid(t)

        rkm, eta, t_mfluid = load_surf_mfluid(outdir_mfluid, frameno_mfluid)
        mfluid_plot.set_data(rkm, eta)
        if max(abs(t_ale-t),abs(t_mfluid-t)) > 0.1:
            print('+++ t = %.1f: frameno_ale = %i at t = %1.f, frameno_mfluid = %i at t = %1.f' \
                % (t,frameno_ale, t_ale, frameno_mfluid, t_mfluid))

        title_text.set_text('Hemispherical crater, radius 3000 m, on 4km ocean,' \
                + '  t = %6.1f sec' % t)


    print('Making anim...')
    times = tf_mfluid[:,1]
    #times = [0, 30, 60]
    anim = animation.FuncAnimation(fig, update, frames=times,
                                   interval=200, blit=False)

    # Output files:
    name = 'RC3000_animation2'

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
