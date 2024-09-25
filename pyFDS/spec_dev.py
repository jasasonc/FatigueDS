import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import scipy
from scipy.special import gamma 
import rainflow
from concurrent.futures import ThreadPoolExecutor

from . import tools

class SpecDev:

    def __init__(self, f0_data, damp=None, Q=None, signal_type=None, b=None, C=None, K=None):
        """
        :param freq_data: tuple containing (f0_start, f0_stop, f0_step) [Hz] or a frequency vector
        :param b: S-N curve slope from Basquin equation
        :param C: material constant from Basquin equation (default: C=1)
        :param K: constant of proportionality between stress and deformation (default: K=1)
        :param damp: damping ratio [/]
        :param Q: damping Q-factor [/]
        """

        self.valid_signal_types = ['sine', 'sine_sweep', 'random']

        #check f0_data input
        if (isinstance(f0_data, (tuple)) and len(f0_data)==3) or (
            isinstance(f0_data, np.ndarray) and f0_data.ndim == 1
        ):
           self.f0_range = tools.get_f0_range(self,f0_data)
        else:
            raise ValueError('`f0` should be a tuple containing (f0_start, f0_stop, f0_step) [Hz] or a frequency vector')
        
        #check damping input (Q or damp)
        if damp is not None or Q is not None:
            tools.convert_Q_damp(self,Q=Q, damp=damp)    
            if damp is not None and Q is not None:
                print('Both `damp` and `Q` are defined. Prioritizing `damp`.')
        else:
            raise ValueError('Either `damp` or `Q` must be provided.')
        
        #set signal type, if given in input
        if signal_type is not None:
            self.set_signal_type(signal_type=signal_type)
        
        if b is not None:
            if isinstance(b, (int, float)):
                self.b = b
            else:
                raise ValueError('`b` should be a scalar value')
        
        if C is not None:
            if isinstance(C, (int, float)):
                self.C = C
            else:
                raise ValueError('`C` should be a scalar value')
        
        if K is not None:
            if isinstance(K, (int, float)):
                self.K = K
            else:
                raise ValueError('`K` should be a scalar value')
        
   
    def set_signal_type(self, signal_type):
        """
        Set the type of the signal and update docstring in get_ers and get_fds methods.

        :param signal_type: signal type

        Available signal types: `sine`, `sine_sweep`, `random`
        """
        if isinstance(signal_type, str) and signal_type in self.valid_signal_types:
                self.signal_type = signal_type
                
                if self.signal_type == 'sine':
                    tools.update_docstring(self.get_ers,self.sine)
                    tools.update_docstring(self.get_fds,self.sine)
                if self.signal_type == 'sine_sweep':
                    tools.update_docstring(self.get_ers,self.sine_sweep)
                    tools.update_docstring(self.get_fds,self.sine_sweep)
                if self.signal_type == 'random':
                    tools.update_docstring(self.get_ers,self.random)
                    tools.update_docstring(self.get_fds,self.random)
        else:
            raise ValueError(f"Invalid signal type. Expected one of the following: {self.valid_signal_types}")
        
    
    def get_ers(self,*args,**kwargs):
        """
        get ERS of a signal

        ---

        docstring of a selected signal type
        """        
        if self.signal_type == 'sine':
            self.ers = self.sine(output='ERS', *args, **kwargs)
        
        if self.signal_type == 'sine_sweep':
            #self.ers = self.sine_sweep(output='ERS', *args, **kwargs)
            pass
        
        if self.signal_type == 'random':
            #self.ers = self.random(output='ERS', *args, **kwargs)        
            pass


    def get_fds(self, *args,**kwargs):
        """
        get FDS of a signal

        ---

        docstring of a selected signal type
        """
        if not any(hasattr(self, attr) for attr in ['b', 'C', 'K']):
            raise ValueError('b, C and K parameters must be provided')
        
        if self.signal_type == 'sine':
            self.fds = self.sine(output='FDS',*args, **kwargs)
        
        if self.signal_type == 'sine_sweep':
            #self.ers = self.sine_sweep(output='FDS', *args, **kwargs)
            pass
        
        if self.signal_type == 'random':
            #self.ers = self.random(output='FDS', *args, **kwargs)        
            pass


    def sine(self, sine_freq=None, amp=None, t_total=None, exc_type='acc', output=None):
        """
        Sine signal

        :param sine_freq: sine frequency [Hz]
        :param amp: signal amplitude [-]
        :param t_total: total excitation time [s] (only for FDS calculation)
        :param exc_type: excitation type (supported: 'acc', 'vel' and 'disp'), default: 'acc'
        :returns: MRS or FDS of defined sine signal
        """
        if all([sine_freq, amp, exc_type]):
            self.sine_freq = sine_freq
            self.amp = amp
            self.exc_type = exc_type

        if t_total is not None:
            self.t_total = t_total

        if not (hasattr(self, 'sine_freq') and hasattr(self, 'sine_freq')):
            raise ValueError('Missing parameter(s). `sine_freq` and `amp` must be provided')

        ω_0i = 2*np.pi*self.f0_range
        
        
        if self.exc_type in ['acc','vel','disp']:
            
            if self.exc_type=='acc':
                a = 0
            elif self.exc_type=='vel':
                a = 1
            elif self.exc_type=='disp':
                a = 2
        else:
            raise ValueError(f"Invalid excitation type. Supported types: `acc`, `vel` and `disp`.")

        #ERS

        if output == 'ERS':

            R_i = -self.amp*(ω_0i)**a/(np.sqrt((1-(self.sine_freq/self.f0_range)**2)**2 + (self.sine_freq/(self.Q*self.f0_range))**2))
            return np.abs(R_i) 

        #FDS
        elif output == 'FDS':

            if not hasattr(self, 't_total'):
                raise ValueError('Missing parameter `t_total`.')
            
            h = self.sine_freq / self.f0_range
            D_i = self.K**self.b/self.C * self.f0_range * self.t_total * self.amp**self.b * ω_0i**(self.b*(a-2)) * h**(a*self.b+1) / ((1-h**2)**2 + (h/self.Q)**2)**(self.b/2)
            return D_i 


    def sine_sweep(self, const_amp, const_f_range, exc_type='acc', dt=1, sweep_type=None, sweep_rate=None):
        pass


    def random():
        pass