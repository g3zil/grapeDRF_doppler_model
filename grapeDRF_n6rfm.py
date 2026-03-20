#!/usr/bin/env python
import os
import datetime
import logging
logger = logging.getLogger(__name__)

import pickle
import numpy as np
from scipy import signal
import matplotlib as mpl
from matplotlib import pyplot as plt
import digital_rf as drf
from eclipse_calc import solarContext

# --- Matplotlib global settings ---
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.linestyle'] = ':'
mpl.rcParams['figure.figsize'] = np.array([15, 8])
mpl.rcParams['axes.xmargin'] = 0

def load_grape_drf(sDate, eDate, data_dir, channel='ch0'):
    """
    Load GrapeDRF data from Digital RF.
    Handles both single- and multi-frequency data without IndexErrors.
    Returns a dictionary containing:
        - bigarray_dct: dict of {frequency: complex samples}
        - latest_meta: metadata for the last sample
        - properties: DRF channel properties
        - timevec_utc: list of datetime objects for each sample
    """
    meta_dir = os.path.join(data_dir, channel, 'metadata')
    do = drf.DigitalRFReader(data_dir)
    dmr = drf.DigitalMetadataReader(meta_dir)

    # Get bounds for DRF data and metadata
    s0, s1 = do.get_bounds(channel)
    first_sample, last_sample = dmr.get_bounds()
    print("metadata bounds are %i to %i" % (first_sample, last_sample))

    start_idx = int(np.uint64(first_sample))
    print('computed start_idx =', start_idx)

    fields = dmr.get_fields()
    print("Available fields are <%s>" % (str(fields)))

    # --- Read first frequency, latitude, longitude from metadata ---
    data_dict = dmr.read(start_idx, start_idx + 2, "center_frequencies")
    freqList = list(data_dict.values())[0]
    print("freq =", freqList[0])

    data_dict = dmr.read(start_idx, start_idx + 2, "lat")
    theLatitude = list(data_dict.values())[0]
    print("Latitude:", theLatitude)

    data_dict = dmr.read(start_idx, start_idx + 2, "long")
    theLongitude = list(data_dict.values())[0]
    print("Longitude:", theLongitude)

    # --- Latest metadata for center frequencies ---
    latest_meta = dmr.read_latest()
    latest_inx = list(latest_meta.keys())[0]
    cntr_freqs = [float(f) for f in latest_meta[latest_inx]['center_frequencies']]  # ensure float
    print("Center frequencies:", cntr_freqs)

    properties = do.get_properties(channel)
    fs = properties['samples_per_second']

    # --- Convert datetime to sample indices for reading DRF data ---
    sinx_0 = drf.util.time_to_sample(sDate, fs)
    sinx_1 = drf.util.time_to_sample(eDate, fs)
    blks = do.get_continuous_blocks(sinx_0, sinx_1, channel)
    for sinx, nsamps in blks.items():
        break  # Take first continuous block

    # --- Initialize bigarray dictionary for all frequencies ---
    bigarray_dct = {cfreq: np.zeros(nsamps, dtype=complex) for cfreq in cntr_freqs}

    # --- Read DRF data block ---
    data = do.read_vector(sinx, nsamps, channel)

    # --- Handle single- vs multi-frequency DRF data ---
    for cfreq_inx, (cfreq, bigarray) in enumerate(bigarray_dct.items()):
        print(f"\nWorking on {cfreq} MHz...")
        # *** Change added here: handle single vs multi-frequency DRF data ***
        # Prevents IndexError when Digital RF file has only one frequency channel
        if data.ndim == 1:
            # Single-frequency DRF files return 1D array, assign directly
            bigarray[:] = data
        else:
            # Multi-frequency DRF files return 2D array (samples × freq), assign column
            bigarray[:] = data[:, cfreq_inx]

    # --- Build result dictionary ---
    result = {
        'bigarray_dct': bigarray_dct,
        'latest_meta': latest_meta[latest_inx],
        'properties': properties,
        'timevec_utc': [drf.util.sample_to_datetime(sinx, fs) + drf.util.samples_to_timedelta(x, fs)
                         for x in range(nsamps)]
    }
    return result

class GrapeDRF:
    """
    Wrapper class for GrapeDRF data with plotting capabilities.
    """
    def __init__(self, sDate, eDate, station,
                 output_dir=os.path.join('output', 'grapeDRF')):
        # --- Ensure output directory exists ---
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # --- Build file names for cached data and PNG output ---
        sDate_str = sDate.strftime('%Y%m%d.%H%M')
        eDate_str = eDate.strftime('%Y%m%d.%H%M')

        event_fname = f'{sDate_str}-{eDate_str}_{station}_grapeDRF'
        png_fname = event_fname + '.png'
        png_fpath = os.path.join(output_dir, png_fname)
        ba_fpath = os.path.join(output_dir, event_fname + '.ba.pkl')
        data_dir = os.path.join('data', 'psws_grapeDRF', station)

        # --- Load DRF data or cached pickle file ---
        if not os.path.exists(ba_fpath):
            result = load_grape_drf(sDate, eDate, data_dir)
            with open(ba_fpath, 'wb') as fl:
                pickle.dump(result, fl)
        else:
            print(f'Using cached file {ba_fpath}...')
            with open(ba_fpath, 'rb') as fl:
                result = pickle.load(fl)

        self.result = result
        self.cfreqs = list(result['bigarray_dct'].keys())
        self.fs = result['properties']['samples_per_second']
        self.sDate = sDate
        self.eDate = eDate
        self.data_dir = data_dir
        self.event_fname = event_fname
        self.output_dir = output_dir
        self.png_fpath = png_fpath
        # --- Custom colormap for spectrogram ---
        self.cmap = mpl.colors.LinearSegmentedColormap.from_list(" ", ["black", "darkgreen", "green", "yellow", "red"])
        self.spectrum_timevec = None

    def plot_figure(self, cfreqs=None, png_fpath=None, **kwargs):
        """
        Plot all specified frequencies (or all if None) to a single figure.
        """
        print(f'Now plotting {self.event_fname}...')
        if cfreqs is None:
            cfreqs = self.cfreqs

        ncols = 1
        nrows = len(cfreqs)
        ax_inx = 0
        fig = plt.figure(figsize=(15, 4 * nrows))
        for cfreq in cfreqs:
            print(f'   {cfreq} MHz...')
            ax_inx += 1
            ax = fig.add_subplot(nrows, ncols, ax_inx)
            self.plot_ax(cfreq, ax, **kwargs)

        fig.tight_layout()
        if png_fpath is None:
            png_fpath = self.png_fpath

        fig.savefig(png_fpath, bbox_inches='tight')
        print(png_fpath)

    def plot_ax(self, cfreq, ax, cmap=None, plot_colorbar=False,
                xlim=None, solar_lat=None, solar_lon=None,
                overlaySolarElevation=True, overlayEclipse=False):
        """
        Plot spectrogram for a single frequency on the provided axis.
        Overlays solar elevation and eclipse if requested.
        """
        ylabel = ['Doppler Shift (Hz)']
        ax.set_ylabel('\n'.join(ylabel))
        ax.set_xlabel('UTC')

        bigarray = self.result['bigarray_dct'].get(cfreq)
        if bigarray is None:
            msg = f'ERROR: No data for {cfreq} MHz'
            ax.text(0.5, 0.5, msg, ha='center', va='center', transform=ax.transAxes)
            print(msg)
            return

        # --- Compute spectrogram ---
        f, t_spec, Sxx = signal.spectrogram(bigarray, fs=self.fs, nfft=1024, window='hann', return_onesided=False)

        # --- Generate consistent time vector for spectrogram ---
        if self.spectrum_timevec is None:
            ts0 = min(self.result['timevec_utc']).timestamp()
            ts1 = max(self.result['timevec_utc']).timestamp()
            ts_vec = np.linspace(ts0, ts1, len(t_spec))
            self.spectrum_timevec = [datetime.datetime.utcfromtimestamp(x) for x in ts_vec]

        f = np.fft.fftshift(f).astype('float64')
        Sxx = np.fft.fftshift(Sxx, axes=0)
        Sxx_db = 10 * np.log10(Sxx)

        if cmap is None:
            cmap = self.cmap

        mpbl = ax.pcolormesh(self.spectrum_timevec, f, Sxx_db, cmap=cmap)

        # --- Optional colorbar ---
        if plot_colorbar:
            try:
                fig = ax.figure
                fig.colorbar(mpbl, label='PSD [dB]')
            except Exception as e:
                print("Failed to add colorbar:", e)

        # --- Overlay solar/eclipse information ---
        sts = solarContext.solarTimeseries(self.sDate, self.eDate, solar_lat, solar_lon)
        odct = {'color': 'white', 'lw': 4, 'alpha': 0.75}
        if overlaySolarElevation:
            sts.overlaySolarElevation(ax, **odct)
        if overlayEclipse:
            sts.overlayEclipse(ax, **odct)

        # --- Format x-axis time labels ---
        if xlim is None:
            xlim = (self.sDate, self.eDate)

        xticks = ax.get_xticks()
        xtkls = []
        for xtk in xticks:
            dt = mpl.dates.num2date(xtk)
            xtkl = dt.strftime('%H:%M')
            xtkls.append(xtkl)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtkls)
