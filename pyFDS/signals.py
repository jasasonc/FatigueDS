import numpy as np
import scipy
from scipy import signal
from scipy.special import gamma
from . import tools
from tqdm import tqdm
import rainflow

def sine(self, output=None):
    """
    Internal function for calculating ERS and FDS of a sine signal.
    """

    omega_0i = 2*np.pi*self.f0_range

    #Getting the ERS with self.get_ers()
    if output == 'ERS':

        R_i = -self.amp*(omega_0i)**self.a/(np.sqrt((1-(self.sine_freq/self.f0_range)**2)**2 + (self.sine_freq/(self.Q*self.f0_range))**2))
        return np.abs(R_i) 


    #Getting the FDS with self.get_fds()
    elif output == 'FDS':

        if not hasattr(self, 't_total'):
            raise ValueError('Missing parameter `t_total`.')

        h = self.sine_freq / self.f0_range
        D_i = self.K**self.b/self.C * self.f0_range * self.t_total * self.amp**self.b * omega_0i**(self.b*(self.a-2)) * h**(self.a*self.b+1) / ((1-h**2)**2 + (h/self.Q)**2)**(self.b/2)
        return D_i 


def sine_sweep(self, output=None):
    """
    Internal function for calculating ERS and FDS of a sine sweep signal.
    """
    
    R_i_all = np.zeros((len(self.f0_range), len(self.const_amp)))
    fds = np.zeros(len(self.f0_range))
    ers = np.zeros(len(self.f0_range))
    
    for i in range(len(self.f0_range)):
        omega_0i = 2 * np.pi * self.f0_range[i]
        
        D_i = 0

        for n in range(len(self.const_amp)):
            amp = self.const_amp[n]
            f1 = self.const_f_range[n]
            f2 = self.const_f_range[n + 1]
            h1 = f1 / self.f0_range[i]
            h2 = f2 / self.f0_range[i]

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
            
            
                const = self.K**self.b / self.C * self.f0_range[i] * tb * amp**self.b * omega_0i**(self.b * (self.a - 2))
                integral = scipy.integrate.trapezoid(M_h * h**(self.a * self.b - 1) / ((1 - h**2)**2 + (h / self.Q)**2)**(self.b / 2), x=h)
                fds[i] += const * integral

            elif output == 'ERS':
                if self.f0_range[i] <= f1:
                    Omega_1 = 2*np.pi*f1
                    R_i = Omega_1**self.a * amp / (np.sqrt((1 - h1**2)**2 + (h1 / self.Q)**2))  # page 32/501 eq. [1.22]
                elif self.f0_range[i] >= f2:
                    Omega_2 = 2*np.pi*f2
                    R_i = Omega_2**self.a * amp / (np.sqrt((1 - h2**2)**2 + (h2 / self.Q)**2))  # page 32/501 eq. [1.23]
                else:
                    R_i = omega_0i**self.a * amp * self.Q  # page 31/501 eq. [1.21]
                
                R_i_all[i, n] = R_i
            


        ers[i] = max(R_i_all[i, :])
    
    if output == 'ERS':
        return ers
    elif output == 'FDS':
        return fds
    

def random_psd(self, output=None):  
        """
        Internal function for calculating ERS and FDS of a random signal in frequency domain.
        """
        
        fds = np.zeros(len(self.f0_range))
        ers = np.zeros(len(self.f0_range))        
        D_i = 0      
        
        # constants
        C0 = np.pi/(4*self.damp)
        C_disp = C0 * 1/((2*np.pi)**4 * self.f0_range**3)
        C_vel = C0 * 1/((2*np.pi)**2 * self.f0_range)
        C_acc = C0 * self.f0_range

        # rms sums
        z_rms_2 = tools.rms_sum(f_0=self.f0_range, psd_freq=self.psd_freq, psd_data=self.psd_data, damp=self.damp, motion='rel_disp') * C_disp
        z_rms = np.sqrt(z_rms_2) * self.unit_scale
        
        dz_rms_2 = tools.rms_sum(f_0=self.f0_range, psd_freq=self.psd_freq, psd_data=self.psd_data, damp=self.damp, motion='rel_vel') * C_vel
        dz_rms = np.sqrt(dz_rms_2) * self.unit_scale
        
        if output == 'FDS': #ddz only needed for FDS calculation
            ddz_rms_2 = tools.rms_sum(f_0=self.f0_range, psd_freq=self.psd_freq, psd_data=self.psd_data, damp=self.damp, motion='rel_acc') * C_acc 
            ddz_rms = np.sqrt(np.abs(ddz_rms_2)) * self.unit_scale

        # ERS calculation
        if output == 'ERS':
            n0 = 1/np.pi * dz_rms/z_rms
            ers = (2*np.pi*self.f0_range)**2 * z_rms * np.sqrt(2*np.log(n0*self.T))
            return ers
        
        # FDS calculation (damage according to Vol. 0, page 89/198, equation (A1-93))
        elif output == 'FDS':
            
            np_plus = 1/(2*np.pi) * ddz_rms/dz_rms
            n0 = 1/np.pi * dz_rms/z_rms
            fds = self.K**self.b/self.C * n0 * self.T * (z_rms*np.sqrt(2))**self.b * gamma(1 + self.b/2)
            return fds


def random_time(self, output=None):
    """
    Internal function for calculating ERS and FDS of a sine random signal in time domain.
    """

    if output=='ERS':
        ers = np.zeros(len(self.f0_range))
        for i in tqdm(range(len(self.f0_range))):               
            z = tools.response_relative_displacement(self.time_data*self.unit_scale, self.dt, f_0=self.f0_range[i], damp=self.damp)
            R_i = np.max(z) * (2*np.pi*self.f0_range[i])**2 
            ers[i] = R_i
        return ers
    
    if output=='FDS':
        fds = np.zeros(len(self.f0_range))
        
        for i in tqdm(range(len(self.f0_range))):                    
            z = tools.response_relative_displacement(self.time_data*self.unit_scale, self.dt, f_0=self.f0_range[i], damp=self.damp)
            
            rf = rainflow.count_cycles(z)
            rf = np.asarray(rf)
            cyc_sum = np.sum(rf[:,1]*2 * (rf[:,0]/2)**self.b) # *2 and /2 because rainflow returns cycles and ranges, fds theory is defined for half cycles and amplitudes
            # if hasattr(self, 'T'):
            #     D_i = self.T / self.t_total * self.K**self.b / (self.C) * cyc_sum
            # else:
            D_i = self.K**self.b / (self.C) * cyc_sum
            fds[i] = D_i
        return fds
    
    #TODO: Implement multiprocessing

    #from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
    
    # def process_frequency(f_0):
    #     z = tools.response_relative_displacement(self.time_data*self.unit_scale, self.dt, f_0=f_0, damp=self.damp)
    #     if output == 'ERS':
    #         R_i = np.max(z) * (2*np.pi*f_0)**2 
    #         return R_i
    #     elif output == 'FDS':
    #         rf = rainflow.count_cycles(z)
    #         rf = np.asarray(rf)
    #         cyc_sum = np.sum(rf[:,1] * rf[:,0]**self.b)
    #         if hasattr(self, 'T'):
    #             D_i = self.T / self.t_total * self.K**self.b / (self.C) * cyc_sum
    #         else:
    #             D_i = self.K**self.b / (self.C) * cyc_sum
    #         return D_i
        
    # with ThreadPoolExecutor() as executor:
    #     print('Calculating fatigue damage for each SDOF system...')
    #     result = np.fromiter(tqdm(executor.map(process_frequency, self.f0_range), total=len(self.f0_range)), dtype=float)
        
    #     return result
        