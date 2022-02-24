import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylab import rcParams

bhnum = [1,2,3,4,5,8,9,16,25]
areas = [1600, 1600, 1600, 900, 400, 100, 400]
numpiles = [25, 9, 4, 16, 9, 4, 4]

#soil parameters
cov = 80
sofs = [1, 10, 20]

#swap the above parameters for these ones if yo uwant to lok at a COV of 40
#cov = 40
#sofs = [10]

folder = 'si_results_avediff/' #the folder you want to look at

plotheight = 3

#read in borehole location and fitness information
foldername = 'si_results_phoenix/si_results/' + folder
coords = list()
fitness = list()
for b in range(len(numpiles)):
    coords.append(list())
    fitness.append(list())
    for sof in sofs:
        xcoords = pd.read_fwf(foldername+'FinalEA-Xcoords_piles-'+str(numpiles[b])+'_BArea-'+str(areas[b])+'_sof-'+str(sof)+'_cov-'+str(cov)+'_anis-1.txt').values
        ycoords = pd.read_fwf(foldername+'FinalEA-Ycoords_piles-'+str(numpiles[b])+'_BArea-'+str(areas[b])+'_sof-'+str(sof)+'_cov-'+str(cov)+'_anis-1.txt').values
        fitness[-1].append(pd.read_fwf(foldername+'FinalEA-stats_piles-'+str(numpiles[b])+'_BArea-'+str(areas[b])+'_sof-'+str(sof)+'_cov-'+str(cov)+'_anis-1.txt',skiprows=5).values[:,7])
        coords[-1].append([xcoords,ycoords])
coords = np.array(coords)
fitness = np.array(fitness)


#read in pile information
pl = list()
for b in range(len(numpiles)):
    pl.append(pd.read_fwf(foldername+'pile_locations-'+str(numpiles[b])+'_BArea-'+str(areas[b])+'.txt',header=None).values)


#plot borehole locations as one big plot for easy comparison

#loop through the scales of fluctuation
for sof in range(len(sofs)):

    fig,plots = plt.subplots(nrows=coords.shape[0], ncols=coords.shape[3],figsize=(plotheight*coords.shape[3],plotheight*coords.shape[0]))

    for i in range(coords.shape[3]):
        plots[0, i].set_title(str(bhnum[i]) + ' Boreholes')
        for j in range(coords.shape[0]):
            plots[j,i].plot(pl[j][0,:],pl[j][1,:],'s',color='r')
            plots[j,i].plot(coords[j,sof,0,i,:bhnum[i]],coords[j,sof,1,i,:bhnum[i]],'o',color='b')
            plots[j,i].set_xlim(0,80)
            plots[j,i].set_ylim(0,80)
            plots[j,i].set_xticks(np.linspace(0,80,5))
            plots[j,i].set_yticks(np.linspace(0,80,5))
        plots[-1,i].set_xlabel('Site length (m)')
    plots[0,-1].legend(('Pile','BH'),numpoints=1,ncol=2, labelspacing=0, handlelength=0.1)
    for j in range(coords.shape[0]):
        plots[j,0].set_ylabel(str(numpiles[j])+ ' piles, '+str(areas[j])+r' m$^2$')
    fig.tight_layout()
    fig.savefig(foldername+'Borehole-locations_'+str(sofs[sof])+'-sof_'+str(cov)+'-cov.png',dpi=300)
    #fig.savefig(location+foldernames[f]+'\\Borehole locations.png',dpi=300)
    plt.close()


    rcParams['figure.figsize'] = 3,3


    #save individual location plots
    for i in range(coords.shape[3]):
        for j in range(coords.shape[0]):
            plt.plot(pl[j][0, :], pl[j][1, :], 's', color='r')
            plt.plot(coords[j, sof, 0, i, :bhnum[i]], coords[j, sof, 1, i, :bhnum[i]], 'o', color='b')
            plt.xlabel('Site length (m)')
            plt.xlim(0, 80)
            plt.ylim(0, 80)
            plt.tight_layout()
            plt.savefig(foldername +'locations_'+str(sofs[sof])+'-sof_'+str(cov)+'-cov_'+ str(bhnum[i]) + '-boreholes_' + str(numpiles[j]) + '-piles_' + str(areas[j]) + '-area.png', dpi=300)
            plt.close()

    rcParams['figure.figsize'] = 5,4

    #save individual fitness plots
    for j in range(coords.shape[0]):
        plt.plot(bhnum,fitness[j,sof,:])
        plt.xlabel('No. Boreholes')
        plt.ylabel('Average diff. set. (m/m)')
        plt.tight_layout()
        plt.savefig(foldername +'fitness_'+str(sofs[sof])+'-sof_'+str(cov)+'-cov_'+ str(numpiles[j]) + '-piles_' + str(areas[j]) + '-area.png', dpi=300)
        plt.close()

    plt.style.use('default')


