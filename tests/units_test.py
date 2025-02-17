import os
import sys
import pytest
import numpy as np

my_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, my_path + '/../')

import FatigueDS
from test_data import *


def test_version():
    """ check sdypy_template_project exposes a version attribute """
    assert hasattr(FatigueDS, "__version__")
    assert isinstance(FatigueDS.__version__, str)

class TestUnits:
    """ Testing ers and fds funtions unit scaling. 
    
    - ERS: no scaling applied
    - FDS: scaling is applied.
    """

    def test_units_time(self):
        time_data = np.load('test_data/test_time_history.npy', allow_pickle=True)
        time_history_data = time_data[:,1]
        t = time_data[:,0] 
        dt = t[2] - t[1]

        sd_g = FatigueDS.SpecificationDevelopment(freq_data=(20, 200, 5))  # time history (psd averaging)
        sd_ms2 = FatigueDS.SpecificationDevelopment(freq_data=(20, 200, 5))  # time history (psd averaging)

        sd_g.set_random_load((time_history_data, dt), unit='g')  # (time history, dt)
        sd_ms2.set_random_load((time_history_data * 9.81, dt), unit='ms2')  # (time history, dt)

        sd_g.get_ers()
        sd_ms2.get_ers()

        sd_g.get_fds(b=5, C=1, K=1)
        sd_ms2.get_fds(b=5, C=1, K=1)

        assert np.allclose(sd_g.ers * 9.81, sd_ms2.ers)
        assert np.allclose(sd_g.fds, sd_ms2.fds)

    def test_units_psd(self):
        _psd_data = np.load('test_data/test_psd.npy', allow_pickle=True)
        psd_freq = _psd_data[:,0]
        psd_data = _psd_data[:,1]

        sd_PSD_g = FatigueDS.SpecificationDevelopment(freq_data=(20, 200, 5))
        sd_PSD_ms2 = FatigueDS.SpecificationDevelopment(freq_data=(20, 200, 5))
        sd_PSD_g.set_random_load((psd_data, psd_freq), unit='g', T=133.5711234541)
        sd_PSD_ms2.set_random_load((psd_data * 9.81**2, psd_freq), unit='ms2', T=133.5711234541)
        sd_PSD_g.get_ers()
        sd_PSD_ms2.get_ers()
        sd_PSD_g.get_fds(b=5, C=1, K=1)
        sd_PSD_ms2.get_fds(b=5, C=1, K=1)

        assert np.allclose(sd_PSD_g.ers * 9.81, sd_PSD_ms2.ers)
        assert np.allclose(sd_PSD_g.fds, sd_PSD_ms2.fds)        
