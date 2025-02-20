#!/bin/env python

import os
import datetime
import logging
logger  = logging.getLogger(__name__)

import numpy as np
import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt

import grapeDRF

letters = 'abcdefghijklmnopqrztuvwxyz'

mpl.rcParams['font.size']       = 12
mpl.rcParams['font.weight']     = 'bold'
mpl.rcParams['axes.grid']       = True
mpl.rcParams['axes.titlesize']  = 30
mpl.rcParams['grid.linestyle']  = ':'
mpl.rcParams['figure.figsize']  = np.array([15, 8])
mpl.rcParams['axes.xmargin']    = 0
mpl.rcParams['legend.fontsize'] = 'xx-large'

station_dct = {}
sdct    = station_dct['w2naf'] = {}
sdct['QTH'] = 'Spring Brook, PA'

if __name__ == '__main__':
    output_dir = os.path.join('output','w2naf_2024eclipse')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sDate       = datetime.datetime(2024,4,8)
    eDate       = datetime.datetime(2024,4,9)
    station     = 'w2naf'
    lat         =  41.335116 # W2NAF
    lon         = -75.600692 # W2NAF

    figd = {}
    figd['solar_lat']               = lat
    figd['solar_lon']               = lon
    figd['overlaySolarElevation']   = True
    figd['overlayEclipse']          = True
    figd['xlim']                    = (sDate,eDate)

    # 'center_frequencies': array([ 2.5 ,  3.33,  5.  ,  7.85, 10.  , 14.67, 15.  , 20.  , 25.  ])
    cfreqs          = [20,15,10,5]
    # cfreqs          = [3.33,7.85,14.67]
    plot_list   = []
    plot_list.append('WDgrape')
    # plot_list.append('VLF')
#    plot_list.append('gmag')

    str_sDate   = sDate.strftime('%Y%m%d.%H%M')
    str_eDate   = eDate.strftime('%Y%m%d.%H%M')
    png_ = []
    png_.append(str_sDate)
    png_.append(str_eDate)
    png_.append(station)
    for pll in plot_list:
        png_.append(pll)
        if pll == 'WDgrape':
            png_ = png_ +  ['{!s}'.format(x) for x in cfreqs]

    png_fname   = '_'.join(png_)+'.png'
    png_fpath   = os.path.join(output_dir,png_fname)


    nrows       = len(plot_list)
    if 'WDgrape' in plot_list:
        nrows += len(cfreqs) - 1
    ncols       = 1
    ax_inx      = 0
    axs         = []

    fig         = plt.figure(figsize=(22,nrows*5))
    letter_fdict = {'size':32}
    # Grape Plots ##########################
    if 'WDgrape' in plot_list:
        gDRF                        = grapeDRF.GrapeDRF(sDate,eDate,station)
        g_figd                      = figd.copy()
        for cfreq in cfreqs:
            print('   {!s} MHz...'.format(cfreq))
            ax_inx      += 1
            ax          = fig.add_subplot(nrows,ncols,ax_inx)
            axs.append(ax)
            gDRF.plot_ax(cfreq,ax,**g_figd)
            ax.set_title('({!s})'.format(letters[ax_inx-1]),loc='left',fontdict=letter_fdict)
            ax.set_title('{!s} MHz Receiver'.format(cfreq))

            rt_dop_csv = os.path.join('data','shibaji',f'LoS_doppler_{cfreq}.csv')
            if os.path.exists(rt_dop_csv):
                print(f'FOUND: {rt_dop_csv}')
                rt_dop = pd.read_csv(rt_dop_csv)
                rt_dop['time'] = pd.to_datetime(rt_dop['time'])

                xx = rt_dop['time']
                yy = rt_dop['fd']

                ax.plot(xx,yy,color='r',lw=5,zorder=1000)
            else:
                print(f'NOT FOUND: {rt_dop_csv}')

    # Finalize Figure ######################
    for ax_inx,ax in enumerate(axs):
        ax.set_xlim(sDate,eDate)
        xticks  = ax.get_xticks()
        ax.set_xticks(xticks)
        if ax_inx != len(axs)-1:
            ax.set_xlabel('')
            xtkls = ['']*len(xticks)
        else:
            ax.set_xlabel('UTC')
            xtkls   = []
            for xtk in xticks:
                dt      = mpl.dates.num2date(xtk)
                xtkl    = dt.strftime('%H:%M')
                xtkls.append(xtkl)
        ax.set_xticklabels(xtkls)

    sdct    = station_dct.get(station,{})
    if 'QTH' in sdct:
        stxt = '{!s} ({!s})'.format(station.upper(),sdct['QTH'])
    else:
        stxt = station.upper()

    txt = []
    txt.append(stxt)
    txt.append(sDate.strftime('%d %b %Y'))
    fontdict    = {'size':42,'weight':'bold'}
    fig.text(0.5,1.,'\n'.join(txt),fontdict=fontdict,ha='center',va='bottom')

    fig.tight_layout()
    fig.savefig(png_fpath,bbox_inches='tight')
    print(png_fpath)
