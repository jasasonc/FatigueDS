
import pytest
import numpy as np
import sys
import os

my_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, my_path + '/../')

import FatigueDS

from test_data import *

# Pytest will discover and run all test functions named `test_*` or `*_test`.

def test_version():
    """ check sdypy_template_project exposes a version attribute """
    assert hasattr(FatigueDS, "__version__")
    assert isinstance(FatigueDS.__version__, str)


class TestCore:
    """ Testing core functions """

    def test_sine(self):
        """ Test the sine function """
        sd_sine = FatigueDS.SpecificationDevelopment(freq_data=(0, 2000, 5))
        sd_sine.set_sine_load(sine_freq=500, amp=10, t_total=3600)
        sd_sine.get_ers()
        sd_sine.get_fds(b=5, C=1, K=1)

        assert np.allclose(sd_sine.ers, sine_ers_true)
        assert np.allclose(sd_sine.fds, sine_fds_true)

    def test_sine_sweep(self):
        """ Test the sine sweep function """
        sd_sine_sweep = FatigueDS.SpecificationDevelopment(freq_data=(0, 2000, 5))
        sd_sine_sweep.set_sine_sweep_load(const_amp=[5,10,20], const_f_range=[20,100,500,1000], exc_type='acc', sweep_type='log', sweep_rate=1)
        sd_sine_sweep.get_ers()
        sd_sine_sweep.get_fds(b=5, C=1, K=1)

        assert np.allclose(sd_sine_sweep.ers, sine_sweep_ers_true)
        assert np.allclose(sd_sine_sweep.fds, sine_sweep_fds_true)

    def test_random_psd(self):
        """ Test the random psd function """
        _psd_data = np.load('test_data/test_psd.npy', allow_pickle=True)
        psd_freq = _psd_data[:,0]
        psd_data = _psd_data[:,1]

        sd_PSD = FatigueDS.SpecificationDevelopment(freq_data=(20, 200, 5))
        sd_PSD.set_random_load((psd_data, psd_freq), unit='g', T=133.5711234541)
        sd_PSD.get_ers()
        sd_PSD.get_fds(b=5, C=1, K=1)

        assert np.allclose(sd_PSD.ers, random_psd_ers_true)
        assert np.allclose(sd_PSD.fds, random_psd_fds_true)
    
    def test_random_time_convolution(self):
        """ Test the random time history function with convolution"""
        _time_data = np.load('test_data/test_time_history.npy', allow_pickle=True)
        time_history_data = _time_data[:,1]
        t = _time_data[:,0] 
        dt = t[2] - t[1]

        sd_convolution = FatigueDS.SpecificationDevelopment(freq_data=(20, 200, 5))
        sd_convolution.set_random_load((time_history_data, dt), unit='g')
        sd_convolution.get_ers()
        sd_convolution.get_fds(b=5, C=1, K=1)

        assert np.allclose(sd_convolution.ers, random_time_convolution_ers_true)
        assert np.allclose(sd_convolution.fds, random_time_convolution_fds_true)
    
    def test_random_time_psd_averaging(self):
        """ Test the random time history function with psd averaging"""
        _time_data = np.load('test_data/test_time_history.npy', allow_pickle=True)
        time_history_data = _time_data[:,1]
        t = _time_data[:,0] 
        dt = t[2] - t[1]

        sd_averaging = FatigueDS.SpecificationDevelopment(freq_data=(20, 200, 5))
        sd_averaging.set_random_load((time_history_data, dt), unit='g', method='psd_averaging', bins=10)
        sd_averaging.get_ers()
        sd_averaging.get_fds(b=5, C=1, K=1)

        assert np.allclose(sd_averaging.ers, random_time_averaging_ers_true)
        assert np.allclose(sd_averaging.fds, random_time_averaging_fds_true)


