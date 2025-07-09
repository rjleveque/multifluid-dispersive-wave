"""
Adjust which outdirs to plot from, and possibly also adjust xlim and ylim
command to show the desired scale based on the size of the crater.

Then run this script from IPython and it should loop through all available
plot frames.  Type a number at the prompt to go to a particular frame,
or return for next frame, q to quit.
"""

from pylab import *
import glob,os,sys

h0 = 4000.
outdir40_thinc = None
outdir20_thinc = None
outdir10_thinc = None
outdir40_muscl = None
outdir20_muscl = '_output_20m_muscl'
outdir10_muscl = '_output_10m_muscl'

# after fixing hllc, gives identical results:
#outdir40_muscl_fix = '_output_40m_muscl_fix'

# which ones to plot:
outdirs = [outdir20_muscl, outdir10_muscl]

# Compare to Airy solution:
outdir_airy = '../../airy_solutions/RC300'
#outdir_airy = None

def plot_airy(outdir_airy, t, color='k'):
    #airy_files = glob.glob('%s/airy*.txt' % outdir_airy)
    #airy_files.sort()
    rkm_airy = None
    eta_airy = None
    fname_airy_pattern = '%s/airy*_t%s.txt' \
                    % (outdir_airy,str(int(times[j])).zfill(4))
    fname_airy = glob.glob(fname_airy_pattern)
    if len(fname_airy) == 0:
        print('Could not find %s' % fname_airy_pattern)
    elif len(fname_airy) == 1:
        d_airy = loadtxt(fname_airy[0])
        rkm_airy = d_airy[:,0]/1e3
        eta_airy = d_airy[:,1]
        plot(rkm_airy, eta_airy, color, label='Airy')
    else:
        print('did not expect fname_airy = ',fname_airy)
    return rkm_airy, eta_airy


def load_mfluid(outdir, j):

    try:
        fortt_files = glob.glob('%s/fort.t*' % outdir)
        fortt_files.sort()
        #print(fortt_files)
        times = []
        for f in fortt_files:
            t = float(open(f).readline().split()[0])
            times.append(t)
        fortc = '%s/fort.c%s' % (outdir, str(j).zfill(4))
        d = loadtxt(fortc)
        rkm = d[:,0] / 1e3
        eta = d[:,2] - h0
        print('frame %i: loaded fort.c at time %.1f sec' % (j,times[j]))
    except:
        print('could not load from ',outdir)
        rkm = None
        eta = None

    return rkm, eta


# find times in fort.t files from the first outdir (assuming all the same):
fortt_files = glob.glob('%s/fort.t*' % outdirs[0])
fortt_files.sort()
#print(fortt_files)
times = []
for f in fortt_files:
    t = float(open(f).readline().split()[0])
    #print('%s: %.2f' % (f,t))
    times.append(t)

# Main loop for plotting...

figure(101, figsize=(8,5))
j = -1
clf()
while j < len(times)-1:
    ans = input('hit return, or int to jump to that frame, or q to quit ')
    if ans == 'q':
        break
    try:
        j = int(ans.replace('j','')) # allow <int> or 'j <int>' as in Iplotclaw
        print('jump to frame %i' % j)
    except:
        j = j+1

    clf()

    c = ['b','r']  # colors for each outdir
    for k,outdir in enumerate(outdirs):
        rkm, eta = load_mfluid(outdir, j)
        label = os.path.split(outdir)[-1]
        label = label.replace('_output_','mfluid ')
        if rkm is not None:
            plot(rkm, eta, color=c[k], label=label)


    if outdir_airy is not None:
        rkm_airy, eta_airy = plot_airy(outdir_airy, times[j], color='k')

    title('time %.0f sec' % times[j])
    grid(True)
    xlim(0,5)
    ylim(-100,300)
    legend(framealpha=1)
    xlabel('radial distance (km)')
    ylabel('surface (m)')
    draw()
