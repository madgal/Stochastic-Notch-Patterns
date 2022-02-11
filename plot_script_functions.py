import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time


def get_data(n, sim, sim_points, n_eq, file):
    '''
    opens file and adds all trajectories together
    discard first n_eq steps
    '''
    # average time to go through a file:
    # 4x4: 3 sec
    # 8x8: 10 sec
    # 16x16: 37 sec

    D = np.zeros((sim * (sim_points - n_eq), n))
    N = np.zeros((sim * (sim_points - n_eq), n))

    t = 0
    nsim = 0
    count = 0
    line_count = 0
    with open(file) as f:
        for line in f:
            t = t + 1
            if t > 1:
                new = line.replace('[', '')
                new = new.replace(']', '')
                new = new.replace(',', '')
                new = new.replace('"', ' ')
                data = np.fromstring(new, sep=' ')
                nsim_new = data[-1]
                if nsim_new != nsim:
                    count = 0
                nsim = nsim_new
                if count > n_eq - 1:
                    D[line_count] = data[0:n]  # t-2
                    N[line_count] = data[2 * n:3 * n]
                    line_count = line_count + 1
                count = count + 1
    return N, D

def get_dataMadisVersion(file,tstart,tend,nsim):

	df = pd.read_csv(filepath_or_buffer=file,header=0).values[tstart:tend,:]

	N,D=[],[]
	fail=False
	try:
		ind = np.argwhere(df[:,3]==nsim)[:,0]
		if len(ind)==0:
			fail=True
		sn = df[ind,2]
		sd = df[ind,0]
		for row in range(len(sn)):
			sn[row] = sn[row].replace('[','').replace(']','')
			sd[row] = sd[row].replace('[','').replace(']','')
			sn[row] = np.array(sn[row].split(",")).astype(float)
			sd[row] = np.array(sd[row].split(",")).astype(float)
			N+= list(sn[row])
			D+= list(sd[row])

	except:
		print simv, "not converted to int"
		fail=True
	if fail:
	   try:
		simv = int(simv)
		ind = np.argwhere(df[:,3]==simv)[:,0]
		sn = df[ind,2]
		sd = df[ind,0]
		for row in range(len(sn)):
			sn[row] = sn[row].replace('[','').replace(']','')
			sd[row] = sd[row].replace('[','').replace(']','')
			sn[row] = np.array(sn[row].split(",")).astype(float)
			sd[row] = np.array(sd[row].split(",")).astype(float)
			print sn[row]
			print N
			N+= list(sn[row])
			print  N
			D+= list(sd[row])
	   except:
		print simv, "not converted to int"

	N = np.array(N)
	D = np.array(D)
	return N,D

def pseudo_potential(N, D, bins):
    '''
    construct the pseudopotential U=-log10(P(N,D))
    from N, D vectors
    '''
    # pipeline to prepare shot noise data for plotting #
    # 1) binning on log scale in [0.1, 10^5]
    # 2) compute 2D histogram of (N, D)
    # 3) change count from 0 to 1 in bins without any count
    # (this is necessary to take the log10 for the pseudopotential)
    # 4) define a 2D matrix 'norm' to account for the uneven binnign when normalizing the PDF
    # 5) to avoid uneven counting of the 'artificial' points with count=1,
    # they are not divided by norm[i][j] but by the maximum possible normalization value norm[-1][-1]

    norm = np.zeros((bins.size - 1, bins.size - 1))
    for i in np.arange(0, bins.size - 1, 1):
        for j in np.arange(0, bins.size - 1, 1):
            norm[i][j] = (bins[i + 1] - bins[i]) * (bins[j + 1] - bins[j])

    dataD = np.ndarray.flatten(D)
    dataN = np.ndarray.flatten(N)
    dist, xbin, ybin = np.histogram2d(dataN, dataD, bins=[bins, bins])

    # add ones to allow computation of U
    dist_new = np.zeros((xbin.size-1, ybin.size-1))
    for i in range(xbin.size - 1):
        for j in range(ybin.size - 1):
            if dist[i][j] == 0:
                dist_new[i][j] = 1.
            else:
                dist_new[i][j] = dist[i][j]

    # normalize probability to 1
    dist_new = dist_new / np.float(np.sum(dist_new))

    # normalize to uneven grid to to get PDF
    for i in range(xbin.size - 1):
        for j in range(ybin.size - 1):
            if dist[i][j] == 0:
                dist_new[i][j] = dist_new[i][j] / norm[-1][-1]
            else:
                dist_new[i][j] = dist_new[i][j] / norm[i][j]
    U = -np.log10(dist_new)
    return U


def plot_landscape(ax, notch, delta, U, bins, xint, yint, bar = True):
    '''
    plot pseudopotential U to axis ax
    '''
    ax.set_xscale('log')
    ax.set_yscale('log')
    color_map = plt.cm.get_cmap('Reds')
    reversed_color_map = color_map.reversed()
    plt.contourf(notch, delta, U, 21, cmap=reversed_color_map)
    if bar:
        cbar = plt.colorbar(orientation='vertical')
        cbar.set_label('Pseudopotential ($-log_{10}P(N,D)$)')
    plt.xlim(xint)
    plt.ylim(yint)
    plt.xlabel('Notch (molecules)')
    plt.ylabel('Delta (molecules)')




##############################################################
#
### functions to find minima and thresholds in pseudopotential
#
###

def find_minima(dist, xbins, notch, delta):
    '''
    find the S and R  minima on pseudopotential landscape
    NB this is custom made already knowing that there are only two main minima
    the condition min(notch[i], delta[j])>1 gets rid of adsorbing boundaries around N=0 or D=0
    that can have a high count in some simulations
    '''
    dist_up = np.ones((xbins.size - 1, xbins.size - 1)) * np.amax(dist)
    dist_low = np.ones((xbins.size - 1, xbins.size - 1)) * np.amax(dist)
    for i in range(xbins.size - 1):
        for j in range(xbins.size - 1):
            if i > j and min(notch[i], delta[j])>1.:
                dist_up[i][j] = dist[i][j]
            elif i < j and min(notch[i], delta[j])>1.:
                dist_low[i][j] = dist[i][j]
    i_s, j_s = np.unravel_index(np.argmin(dist_up), dist_up.shape)
    i_r, j_r = np.unravel_index(np.argmin(dist_low), dist_low.shape)
    return i_s, j_s, i_r, j_r


def find_thresholds(x, U, i_s, j_s, i_r, j_r):
    '''
    find distances from S and R minima where the probability is decreased by a 10-fold
    thus the pseudopotential increases by a unit
    '''
    ### find thresholds:
    # from Receiver: notch moves left, delta moves up
    Uref = U[i_s][j_s]
    i = i_s - 1
    while (U[i][j_s] - Uref) < 1.:
        i = i - 1
    notch_thr_R = x[i]
    j = j_s + 1
    while (U[i_s][j] - Uref) < 1.:
        j = j + 1
    delta_thr_R = x[j]
    # from Sender: notch moves right, delta moves down
    Uref = U[i_r][j_r]
    i = i_r + 1
    while (U[i][j_r] - Uref) < 1.:
        i = i + 1
    notch_thr_S = x[i]
    j = j_r - 1
    while (U[i_r][j] - Uref) < 1.:
        j = j - 1
    delta_thr_S = x[j]
    return notch_thr_R, delta_thr_R, notch_thr_S, delta_thr_S


##############################################################
#
### plot landscapes with minima, fixed points and thresholds
#
###

def plot_to_axis(ax, bins, U, title, printout = False, legend=False):
    '''
    plots pseudopotential U to axis ax
    finds and plots landscape minima corresponding to S and R
    plots fixed points S and R in zero-noise limit
    estimate and plot the threshold for 10-fold probability change
    '''

    x = np.zeros(bins.size - 1)
    for i in range(x.size):
        x[i] = (bins[i] + bins[i + 1]) / 2.

    plot_landscape(ax, bins[0:-1], bins[0:-1], np.transpose(U), bins, [10, 30000], [0.1, 10000])

    # find minima
    i_s, j_s, i_r, j_r = find_minima(U, bins, x, x)
    plt.plot([bins[i_s], bins[i_r]], [bins[j_s], bins[j_r]], 'g*', markersize=10, label='S, R minima')
    # plot the coordinates of S and R in the checkerboard without noise
    plt.plot([505, 5385], [1758, 25], 'k*', markersize=10, label='Fixed points\n($\\sigma=0$)')
    plt.title(title)

    # find and plot thresholds
    notch_thr_R, delta_thr_R, notch_thr_S, delta_thr_S = find_thresholds(x, U, i_s, j_s, i_r, j_r)
    plt.plot([notch_thr_R, 30000.], [delta_thr_R, delta_thr_R], 'k--', label='Threshold for\n10-fold probability\nchange')
    plt.plot([notch_thr_R, notch_thr_R], [0.1, delta_thr_R], 'k--')
    plt.plot([10, notch_thr_S], [delta_thr_S, delta_thr_S], 'k--')
    plt.plot([notch_thr_S, notch_thr_S], [delta_thr_S, 10000], 'k--')

    if legend:
        plt.legend(loc='lower left')

    if printout:
        print(bins[i_r], notch_thr_S, notch_thr_R, bins[i_s])
        print(bins[j_r], delta_thr_S, delta_thr_R, bins[j_s])

