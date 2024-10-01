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

    def __init__(self, freq_data = (10,2000,5), damp=None, Q=10, signal_type=None):
        """
        :param freq_data: tuple containing (f0_start, f0_stop, f0_step) [Hz] or a frequency vector
        :param b: S-N curve slope from Basquin equation
        :param C: material constant from Basquin equation (default: C=1)
        :param K: constant of proportionality between stress and deformation (default: K=1)
        :param damp: damping ratio [/]  
        :param Q: damping Q-factor [/]
        """

        self.valid_signal_types = ['sine', 'sine_sweep', 'random']

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
        
        #set signal type, if given in input
        if signal_type is not None:
            self.set_signal_type(signal_type=signal_type)
        
        # if b is not None:
        #     if isinstance(b, (int, float)):
        #         self.b = b
        #     else:
        #         raise ValueError('`b` should be a scalar value')
        
        # if C is not None:
        #     if isinstance(C, (int, float)):
        #         self.C = C
        #     else:
        #         raise ValueError('`C` should be a scalar value')
        
        # if K is not None:
        #     if isinstance(K, (int, float)):
        #         self.K = K
        #     else:
        #         raise ValueError('`K` should be a scalar value')
        
   
    def set_signal_type(self, signal_type):
        """
        Set the type of the signal and update docstring in get_ers and get_fds methods.

        :param signal_type: Signal type. Available signal types: `sine`, `sine_sweep`, `random`
        """
        if isinstance(signal_type, str) and signal_type in self.valid_signal_types:
                self.signal_type = signal_type
                
                if self.signal_type == 'sine':
                    tools.update_docstring(self.set_load,self.sine)
                if self.signal_type == 'sine_sweep':
                    tools.update_docstring(self.set_load,self.sine_sweep)
                if self.signal_type == 'random':
                    tools.update_docstring(self.set_load,self.random)
        else:
            raise ValueError(f"Invalid signal type. Expected one of the following: {self.valid_signal_types}")
        

    def set_load(self, *args, **kwargs):
        """
        Set the load parameters for the signal type.

        --- 

        docstring of a selected signal type
        """
        if self.signal_type == 'sine':
            self.sine(*args, **kwargs)
        
        if self.signal_type == 'sine_sweep':
            self.sine_sweep(*args, **kwargs)
            
        
        if self.signal_type == 'random':
            #self.random(*args, **kwargs)
            pass



    def get_ers(self,*args,**kwargs):
        """
        get ERS of a signal

        ---

        docstring of a selected signal type
        """        
        if self.signal_type == 'sine':
            self.ers = self.sine(output='ERS', *args, **kwargs)
        
        if self.signal_type == 'sine_sweep':
            self.ers = self.sine_sweep(output='ERS', *args, **kwargs)
            pass
        
        if self.signal_type == 'random':
            #self.ers = self.random(output='ERS', *args, **kwargs)        
            pass


    def get_fds(self, b, C, K, t_total, *args,**kwargs):
        """
        get FDS of a signal

        ---

        docstring of a selected signal type
        """
        if all(isinstance(attr, (int, float)) for attr in [b, C, K]):
            self.b = b
            self.C = C
            self.K = K
        else:
            raise ValueError('b, C and K parameters must be provided')
        
        if isinstance(t_total, (int, float)):
            self.t_total = t_total
        else:
            raise ValueError('missing parameter `t_total`')

        if self.signal_type == 'sine':
            self.fds = self.sine(output='FDS',*args, **kwargs)
        
        if self.signal_type == 'sine_sweep':
            self.fds = self.sine_sweep(output='FDS', *args, **kwargs)
            pass
        
        if self.signal_type == 'random':
            #self.fds = self.random(output='FDS', *args, **kwargs)        
            pass




    def sine(self, sine_freq=None, amp=None, exc_type='acc', output=None):
        """
        Sine signal

        :param sine_freq: sine frequency [Hz]
        :param amp: signal amplitude [m/s^2, m/s, m]
        :param exc_type: excitation type (supported: 'acc [m/s^2]', 'vel[m/s]' and 'disp[m]')
        """

        #Setting the load parameters with self.set_load()
        if output is None:
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


        #Getting the ERS with self.get_ers()
        if output == 'ERS':

            omega_0i = 2*np.pi*self.f0_range

            R_i = -self.amp*(omega_0i)**self.a/(np.sqrt((1-(self.sine_freq/self.f0_range)**2)**2 + (self.sine_freq/(self.Q*self.f0_range))**2))
            return np.abs(R_i) 


        #Getting the FDS with self.get_fds()
        elif output == 'FDS':

            if not hasattr(self, 't_total'):
                raise ValueError('Missing parameter `t_total`.')
            
            omega_0i = 2*np.pi*self.f0_range

            h = self.sine_freq / self.f0_range
            D_i = self.K**self.b/self.C * self.f0_range * self.t_total * self.amp**self.b * omega_0i**(self.b*(self.a-2)) * h**(self.a*self.b+1) / ((1-h**2)**2 + (h/self.Q)**2)**(self.b/2)
            return D_i 


    def sine_sweep(self, const_amp=None, const_f_range=None, exc_type='acc', dt=1, sweep_type=None, sweep_rate=None, output=None):
        """
        Sine sweep signal
        
        :param const_amp: constant amplitude ranges  [m/s^2, m/s, m]
        :param const_f_range: constant frequency ranges [Hz]
        :param exc_type: excitation type (supported: 'acc [m/s^2]', 'vel[m/s]' and 'disp[m]')
        :param dt: time step [s] (default 1 second)
        :param sweep_type: sine sweep type (['linear','lin'] or ['logarithmic','log']) 
        :param sweep_rate: sinusoidal sweep rate [Hz/min] for 'linear' and [oct./min] for 'logarithmic' sweep type
        """

        #Setting the load parameters with self.set_load()
        if output is None:
             if all([const_amp, const_f_range, exc_type, dt, sweep_type, sweep_rate]):
                #neccesarly parameters
                self.const_amp = const_amp
                self.const_f_range = const_f_range
                self.sweep_type = sweep_type
                self.sweep_rate = sweep_rate
                #optional parameters
                self.exc_type = exc_type
                self.dt = dt

        if self.exc_type in ['acc','vel','disp']:
                
            if self.exc_type=='acc':
                self.a = 0
            elif self.exc_type=='vel':
                self.a = 1
            elif self.exc_type=='disp':
                self.a = 2
        else:
            raise ValueError(f"Invalid excitation type. Supported types: `acc`, `vel` and `disp`.")


        #gettng the ERS and FDS with self.get_ers() and self.get_fds()
        f0_len = len(self.f0_range)
        amp_len = len(self.const_amp)
        
        R_i_all = np.zeros((f0_len, amp_len))
        fds = np.zeros(f0_len)
        ers = np.zeros(f0_len)
        
        for i in range(f0_len):
            omega_0i = 2 * np.pi * self.f0_range[i]
            
            D_i = 0
            for n in range(amp_len):
                amp = self.const_amp[n]
                f1 = self.const_f_range[n]
                f2 = self.const_f_range[n + 1]
                h1 = f1 / self.f0_range[i]
                h2 = f2 / self.f0_range[i]
                
                if self.f0_range[i] <= f1:
                    Omega_1 = 2 * np.pi * f1
                    h1 = f1 / self.f0_range[i]
                    R_i = Omega_1**self.a * amp / (np.sqrt((1 - h1**2)**2 + (h1 / self.Q)**2))  # page 32/501 eq. [1.22]
                elif self.f0_range[i] >= f2:
                    Omega_2 = 2 * np.pi * f2
                    h2 = f2 / self.f0_range[i]
                    R_i = Omega_2**self.a * amp / (np.sqrt((1 - (h2)**2)**2 + (h2 / self.Q)**2))  # page 32/501 eq. [1.23]
                else:
                    omega_0i = 2 * np.pi * self.f0_range[i]
                    R_i = omega_0i**self.a * amp * self.Q  # page 31/501 eq. [1.21]
                
                R_i_all[i, n] = R_i
                
                if output == 'FDS':
                    if self.sweep_type is None:
                        raise ValueError("You need to provide either ['linear','lin'] or ['logarithmic','log'] sweep_type.")
                    elif self.sweep_type in ['lin', 'linear']:
                        tb = (self.const_f_range[-1] - self.const_f_range[0]) / self.sweep_rate * 60  # sinusoidal sweep time [s] -> from [Hz/min]
                        dh = (f2 - f1) * self.dt / (self.f0_range[i] * tb)
                        h = np.arange(h1, h2, dh)
                        M_h = h**2 / (h2 - h1)
                    elif self.sweep_type in ['log', 'logarithmic']:
                        tb = 60 * np.log(self.const_f_range[-1] / self.const_f_range[0]) / (self.sweep_rate * np.log(2))  # logarithmic sweep time [s] -> from [oct./min]
                        t = np.arange(0, tb, self.dt)
                        T1 = tb / np.log(h2 / h1)
                        f_t = f1 * np.exp(t / T1)
                        dh = f1 / (T1 * self.f0_range[i]) * np.exp(t / T1) * self.dt
                        h = f_t / self.f0_range[i]
                        M_h = h / (np.log(h2 / h1))
                    else:
                        raise ValueError(f"Invalid method `method`='{self.sweep_type}'. Supported sweep types: 'lin' and 'log'.")
                
                
                    const = self.K**self.b / self.C * self.f0_range[i] * tb * self.const_amp[n]**self.b * omega_0i**(self.b * (self.a - 2))
                    integral = scipy.integrate.trapezoid(M_h * h**(self.a * self.b - 1) / ((1 - h**2)**2 + (h / self.Q)**2)**(self.b / 2), x=h)
                    D_i += const * integral
                
                    fds[i] = D_i

            ers[i] = max(R_i_all[i, :])
        
        if output == 'ERS':
            return ers
        elif output == 'FDS':
            return fds
        
