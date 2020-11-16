import matplotlib.pyplot as plt
import pandas as pd
import sys
from matplotlib import cm
from matplotlib import colors
import colorsys

def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

cmap = cm.get_cmap('hsv')
rgba = cmap(0.5)
all_results = pd.read_csv(sys.argv[1])
NN = all_results["N"]
lightcolors = [ scale_lightness(x,1.75) for x in [colors.ColorConverter.to_rgb(cmap(0.7)),colors.ColorConverter.to_rgb(cmap(0.8)),colors.ColorConverter.to_rgb(cmap(0.9))]]
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(4.8,7))
l1 = ax1.plot(NN,all_results["mean_tau_0.5"],color=cmap(0.7),label=r'0.5')[0]
ax1.fill_between(NN,all_results["low_tau_0.5"],all_results["high_tau_0.5"],color=lightcolors[0])
l2 = ax1.plot(NN,all_results["mean_tau_1.0"],color=cmap(0.8))[0]
ax1.fill_between(NN,all_results["low_tau_1.0"],all_results["high_tau_1.0"],color=lightcolors[1])
l3 = ax1.plot(NN,all_results["mean_tau_1.5"],color=cmap(0.9))[0]
ax1.fill_between(NN,all_results["low_tau_1.5"],all_results["high_tau_1.5"],color=lightcolors[2])
ax1.set_xlabel('N')
ax1.set_ylabel(r'Kendall $\tau$')
ax2.plot(NN,all_results["mean_rmse_0.5"],color=cmap(0.7))
#plt.plot(NN,all_results[:,1],color='black')
#plt.plot(NN,all_results[:,2],color='black')
ax2.fill_between(NN,all_results["low_rmse_0.5"],all_results["high_rmse_0.5"],color=lightcolors[0])
ax2.plot(NN,all_results["mean_rmse_1.0"],color=cmap(0.8))
ax2.fill_between(NN,all_results["low_rmse_1.0"],all_results["high_rmse_1.0"],color=lightcolors[1])
ax2.plot(NN,all_results["mean_rmse_1.5"],color=cmap(0.9))
ax2.fill_between(NN,all_results["low_rmse_1.5"],all_results["high_rmse_1.5"],color=lightcolors[2])
ax2.set_xlabel('N')
ax2.set_ylabel('RMSE [kcal/mol]')
#plt.text(0.1,1.0,r"$\sigma = 1.0$",rotation=45,verticalalignment='center')
line_labels = [r'$\sigma=0.5$',r'$\sigma=1.0$',r'$\sigma=1.5$']
fig.legend([l1, l2, l3],     # The line objects
           labels=line_labels,   # The labels for each line
           loc=(0.6,0.6),   # Position of legend
           borderaxespad=0.1,    # Small spacing around legend box
           #title="Legend Title"  # Title for the legend
           )
#plt.subplots_adjust(hspace=0.25)
plt.savefig('N-CI.pdf',bbox_inches='tight')
