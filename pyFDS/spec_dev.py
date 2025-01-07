import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import scipy
from scipy.special import gamma 
import rainflow
from concurrent.futures import ThreadPoolExecutor

from . import tools
from . import signals

class SpecificationDevelopment:

    def __init__(self, freq_data = (10,2000,5), damp=None, Q=10):
        """
        :param freq_data: tuple containing (f0_start, f0_stop, f0_step) [Hz] or a frequency vector
        :param b: S-N curve slope from Basquin equation
        :param C: material constant from Basquin equation (default: C=1)
        :param K: constant of proportionality between stress and deformation (default: K=1)
        :param damp: damping ratio [/]  
        :param Q: damping Q-factor [/]
        """

        #check freq_data input
        if (isinstance(freq_data, (tuple)) and len(freq_data)==3) or (
            isinstance(freq_data, np.ndarray) and freq_data.ndim == 1
        ):
           self.f0_range = tools.get_freq_range(self,freq_data)
        else:
            raise ValueError('`f0` should be a tuple containing (f0_start, f0_stop, f0_step) [Hz] or a frequency vector')
        
        #check damping input (Q or damp)
        if isinstance(damp, (int, float)) or isinstance(Q, (int, float)):
            tools.convert_Q_damp(self,Q=Q, damp=damp)


    def set_sine_load(self, sine_freq=None, amp=None, exc_type='acc'):
        """
        Sine signal

        :param sine_freq: sine frequency [Hz]
        :param amp: signal amplitude [m/s^2, m/s, m]
        :param exc_type: excitation type (supported: 'acc [m/s^2]', 'vel[m/s]' and 'disp[m]')
        """

        self.signal_type = 'sine'

        if all([sine_freq, amp, exc_type]):
            self.sine_freq = sine_freq
            self.amp = amp
            self.exc_type = exc_type
        else:    
            raise ValueError('Missing parameter(s). `sine_freq` and `amp` must be provided')
            
        if self.exc_type in ['acc','vel','disp']:            
            if self.exc_type=='acc':
                self.a = 0
            elif self.exc_type=='vel':
                self.a = 1
            elif self.exc_type=='disp':
                self.a = 2
        else:
            raise ValueError(f"Invalid excitation type. Supported types: `acc`, `vel` and `disp`.")



    def set_sine_sweep_load(self, const_amp=None, const_f_range=None, exc_type='acc', dt=1, sweep_type=None, sweep_rate=None, ):
        """
        Sine sweep signal
        
        :param const_amp: constant amplitude ranges  [m/s^2, m/s, m]
        :param const_f_range: constant frequency ranges [Hz]
        :param exc_type: excitation type (supported: 'acc [m/s^2]', 'vel[m/s]' and 'disp[m]')
        :param dt: time step [s] (default 1 second)
        :param sweep_type: sine sweep type (['linear','lin'] or ['logarithmic','log']) 
        :param sweep_rate: sinusoidal sweep rate [Hz/min] for 'linear' and [oct./min] for 'logarithmic' sweep type
        """
        
        self.signal_type = 'sine_sweep'
        if all([const_amp, const_f_range, exc_type, dt, sweep_type, sweep_rate]):
            #neccesarly parameters
            self.const_amp = const_amp
            self.const_f_range = const_f_range
            self.sweep_type = sweep_type
            self.sweep_rate = sweep_rate
            #optional parameters
            self.exc_type = exc_type
            self.dt = dt
        else:
            raise ValueError('Missing parameter(s). `const_amp`, `const_f_range`, `sweep_type` and `sweep_rate` must be provided')

        if self.exc_type in ['acc','vel','disp']:   
            if self.exc_type=='acc':
                self.a = 0
            elif self.exc_type=='vel':
                self.a = 1
            elif self.exc_type=='disp':
                self.a = 2
        else:
            raise ValueError(f"Invalid excitation type. Supported types: `acc`, `vel` and `disp`.")  
        

    def set_random_load(self, signal_data=None, T=None, unit='g', method='convolution', bins = None):
        """
        Random signal
        """

        # Signal data must be a tuple
        if isinstance(signal_data, tuple) and len(signal_data) == 2:
        
        # If input is time signal
            if isinstance(signal_data[0], np.ndarray) and isinstance(signal_data[1], (int, float)):
                self.signal_type = 'random_time'
                self.time_data = signal_data[0]  # time-history
                self.dt = signal_data[1] # Sampling interval

                if method in ['convolution', 'psd_averaging']:
                    self.method = method

                else:
                    raise ValueError('Invalid method. Supported methods: `convolution` and `psd_averaging`')
                
                if isinstance(bins, int):
                    self.bins = bins
                if isinstance(T, (int, float)):
                    self.T = T
                    #print('Time duration `T` is not needed for random time signal')
        
        # If input is PSD
            elif isinstance(signal_data[0], np.ndarray) and isinstance(signal_data[1], np.ndarray):
                self.signal_type = 'random_psd'
                self.psd_data = signal_data[0]
                self.psd_freq = signal_data[1]
                
                if isinstance(T, (int, float)):
                    self.T = T
                else:
                    raise ValueError('Time duration `T` must be provided')

            else:
                raise ValueError('Invalid input. Expected a tuple containing (time history data, fs) or (psd data, frequency vector)')
            

        
        # units (to return MRS in m/s^2)
        if unit=='g':
            self.unit_scale = 9.81
        elif unit=='ms2':
            self.unit_scale = 1
        else:
            raise ValueError("Invalid unit selected. Supported units: 'g' and 'ms2'.")


    def get_ers(self):
        """
        get ERS of a signal

        """        
        if self.signal_type == 'sine':
            self.ers = signals.sine(self,output='ERS')
        
        if self.signal_type == 'sine_sweep':
            self.ers = signals.sine_sweep(self,output='ERS')
        
        if self.signal_type == 'random_psd':
            self.ers = signals.random_psd(self,output='ERS')
        
        if self.signal_type == 'random_time':
            if self.method == 'convolution':
                self.ers = signals.random_time(self,output='ERS')   
            elif self.method == 'psd_averaging':
                tools.psd_averaging(self)
                self.ers = signals.random_psd(self,output='ERS')
                


    def get_fds(self,  b, C=1, K=1, t_total=None):
        """
        get FDS of a signal

        :param b: S-N curve slope from Basquin equation
        :param C: material constant from Basquin equation (default: C=1)
        :param K: constant of proportionality between stress and deformation (default: K=1)

        """
        try:
            delattr(self, 't_total')
        except:
            pass

        if isinstance(t_total, (int, float)):
            self.t_total = t_total
        else:
            pass #parameter not needed for all signal types, exception handled in signals.py
        
        if all(isinstance(attr, (int, float)) for attr in [b, C, K]):
            self.b = b
            self.C = C
            self.K = K
        else:
            raise ValueError('Material parameters: b, C and K must be provided')
        
        if self.signal_type == 'sine':
            self.fds = signals.sine(self,output='FDS')
        
        if self.signal_type == 'sine_sweep':
            self.fds = signals.sine_sweep(self,output='FDS')
        
        if self.signal_type == 'random_psd':
            self.fds = signals.random_psd(self,output='FDS')

        if self.signal_type == 'random_time':
            if self.method == 'convolution':
                self.fds = signals.random_time(self,output='FDS')   
            elif self.method == 'psd_averaging':
                tools.psd_averaging(self)
                self.fds = signals.random_psd(self,output='FDS')





        
