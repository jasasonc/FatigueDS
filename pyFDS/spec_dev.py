import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.special import gamma 
import rainflow

from . import tools
from . import signals

class SpecificationDevelopment:

    def __init__(self, freq_data = (10,2000,5), damp=None, Q=10):
        """
        Initialize the SpecificationDevelopment class. Frequency range and damping ratio/Q-factor must be provided.
        Only one of the damping ratio or Q-factor must be provided. If both are provided, damping ratio will be used. If None, Q=10 will be used.

        :param freq_data: tuple containing (f0_start, f0_stop, f0_step) [Hz] or a frequency vector, defining the range where the ERS and FDS will be calculated
        :param damp: damping ratio [/]
        :param Q: damping Q-factor [/] (default: Q=10)
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


    def set_sine_load(self, sine_freq=None, amp=None, t_total=None, exc_type='acc'):
        """
        Set sine signal load parameters

        :param sine_freq: sine frequency [Hz]
        :param amp: signal amplitude [m/s^2, m/s, m]
        :param t_total: total time duration of the signal [s] (only needed for fds calculation)
        :param exc_type: excitation type (supported: 'acc [m/s^2]', 'vel[m/s]' and 'disp[m]')
        """

        self.signal_type = 'sine'

        if all([sine_freq, amp, exc_type]):
            self.sine_freq = sine_freq
            self.amp = amp
            self.exc_type = exc_type
        else:    
            raise ValueError('Missing parameter(s). `sine_freq` and `amp` must be provided')
        
        if isinstance(t_total, (int, float)):
            self.t_total = t_total

            
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
        Set sine sweep signal load parameters
        
        :param const_amp: constant amplitude ranges  [m/s^2, m/s, m]
        :param const_f_range: constant frequency ranges [Hz]
        :param exc_type: excitation type (supported: 'acc [m/s^2]', 'vel[m/s]' and 'disp[m]')
        :param dt: time step [s] (default 1 second)
        :param sweep_type: sine sweep type (['linear','lin'] or ['logarithmic','log']) 
        :param sweep_rate: sinusoidal sweep rate [Hz/min] for 'linear' and [oct./min] for 'logarithmic' sweep type
        """
        
        self.signal_type = 'sine_sweep'
        if all([const_amp, const_f_range, exc_type, dt, sweep_type, sweep_rate]):
            #necessary parameters
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
        

    def set_random_load(self, signal_data=None, T=None, unit='ms2', method='convolution', bins = None):
        """
        Set random signal load parameters

        :param signal_data: tuple containing (time history data, dt) or (psd data, frequency vector)
        :param T: time duration [s]
        :param unit: unit of the signal (supported: 'g' and 'ms2') Parameter only needed for fds calculation
        :param method: method to calculate ERS and FDS (supported: 'convolution' and 'psd_averaging'). Only needed for random time signal
        :param bins: number of bins for PSD averaging method. Only neede for psd averaging method
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
                    print('Time duration `T` is not needed for random time signal')
                self.T = len(self.time_data)*self.dt
        
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
            

        if unit=='g':
            self.unit_scale = 9.81
        elif unit=='ms2':
            self.unit_scale = 1
        else:
            raise ValueError("Invalid unit selected. Supported units: 'g' and 'ms2'.")


    def get_ers(self):
        """
        get extreme response spectrum (ERS) of a signal.

        The unit of the ERS corresponds to the unit of the signal, no scaling is applied.

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
                


    def get_fds(self,  b, C=1, K=1):
        """
        get fatigue damage spectrum (FDS) of a signal.

        Correct unit must be selected in set_random_load method. If unit is 'g', signal is scaled to m/s^2 before FDS calculation, because the FDS theory is based on SI base units.
 
        :param b: S-N curve slope from Basquin equation
        :param C: material constant from Basquin equation (default: C=1)
        :param K: constant of proportionality between stress and deformation (default: K=1)
        """
        
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





        
