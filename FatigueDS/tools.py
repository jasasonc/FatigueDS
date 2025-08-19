import numpy as np
from scipy import signal
from FLife.tools import basquin_to_sn

def convert_Q_damp(self, Q=None, damp=None): 
    """
    Function for converting damping ratio to Q-factor and vice versa.

    :param Q: damping Q-factor [/]
    :param damp: damping ratio [/]
    """

    if damp is not None:
        self.damp = damp
        self.Q = 1 / (2 * self.damp)

    elif Q is not None:
        self.Q = Q
        self.damp = 1 / (2 * self.Q)

def get_freq_range(self, freq_data):
    """
    Function for generating frequency ranges-> X-axis of MRS/FDS plot from freq_data tuple.

    :param freq_data: frequency data (tuple)

    :return: frequency range
    """
    if isinstance(freq_data, tuple) and len(freq_data) == 3:
        f0_start, f0_stop, f0_step = freq_data
        f0_range = np.arange(f0_start, f0_stop + f0_step, f0_step, dtype=float)
    else:
        f0_range = freq_data       
    
    if f0_range[0] == 0:
        f0_range[0] = 1e-3    # sets frequency to a small number to avoid dividing by 0
    
    return f0_range


def rms_sum(f_0, psd_freq, psd_data, damp, motion='rel_disp'):
    """
    This function calculates the response RMS (either relative displacement, velocity or acceleration) for a given 
    natural frequency and damping ratio. 

    :param f_0: system natural frequency [Hz]
    :param psd_freq: PSD frequency range [Hz]
    :param psd_data: PSD data [(m/s^2)^2/Hz] or [g^2/Hz]
    :param damp: damping ratio [/]
    :param motion: which rms sum to perform (supported: rel_disp, rel_vel and rel_acc)

    :return: RMS sum value
    """
    
    df = np.diff(psd_freq)[0]
    rms_sum = 0

    f1 = psd_freq - df / 2
    f2 = psd_freq + df / 2

    # Adjust first and last elements for f1 and f2 respectively
    f1[0] = psd_freq[0]
    f2[-1] = psd_freq[-1]

    for j in range(len(psd_data)):

        h1 = f1[j] / f_0
        h2 = f2[j] / f_0

        # Case where the excitation is defined by PSD comprising "n" straight line segments (Vol.3, equation [8.86])
        
        
        

        if motion == 'rel_disp':
            z_rms = psd_data[j] * (integrals_b(h=h2, b=0, damp=damp) - integrals_b(h=h1, b=0, damp=damp))
            rms_sum += z_rms

        elif motion == 'rel_vel':
            dz_rms = psd_data[j] * (integrals_b(h=h2, b=2, damp=damp) - integrals_b(h=h1, b=2, damp=damp))
            rms_sum += dz_rms
            
        elif motion == 'rel_acc':
            ddz_rms = psd_data[j] * (integrals_b(h=h2, b=4, damp=damp) - integrals_b(h=h1, b=4, damp=damp))
            rms_sum += ddz_rms
       
    return rms_sum



def integrals_b(h, b, damp):
    """
    This function calculates integrals I_b described in [3] and [4]. See equations (A1-74), (A1-75), (A1-76) in [3]
    or [A6.20], [A6.22], [A6.24] in [4] or [8.52], [8.53], [8.54] [4].

    Literature:
        [3] Mechanical Environment Test Specification Development Method - Christian LALANNE
        [4] Christian Lalanne(auth.) Random Vibration Mechanical Vibration and Shock Analysis, Volume 3, Second Edition
    
    :param h: frequency ratio (frequency vs natural frequency) [/]
    :param b: exponent b [/]
    :param damp: damping ratio [/]

    :return: I_b integral value
    """
    
    # constants
    alpha = 2 * np.sqrt(1 - damp**2)    
    beta = 2 * (1 - 2 * damp**2)
    
    C0 = damp / (np.pi * alpha)
    C1 = (h**2 + alpha * h + 1)/(h**2 - alpha * h + 1)
    C2 = (2 * h + alpha) / (2 * damp)
    C3 = (2 * h - alpha) / (2 * damp)
    C4 = 4 * damp / np.pi
    C5 = np.arctan(C2) + np.arctan(C3)
    
    # integrals
    if b == 0:
        Ib = C0 * np.log(C1) + 1 / np.pi * C5  # 84/198 eq. (A1-74) and 560/610 eq. [A6.20]
    
    elif b == 2:
        Ib = C0*np.log(1 / C1) + 1 / np.pi * C5  # 84/198 eq. (A1-75) and 560/610 eq. [A6.22] 
    
    elif b == 4:
        I0 = C0 * np.log(C1) + 1 / np.pi * C5 
        I2 = -C0 * np.log(C1) + 1 / np.pi * C5
        Ib = C4 * h + beta * I2 - I0  # 84/198 eq. (A1-76) and 560/610 eq. [A6.24]

    else:
        raise ValueError(f"Invalid exponent ``b``='{b}'. Supported exponents: 0, 2 and 4.")
    
    return Ib


def response_relative_displacement(time_data, dt, f_0, damp):
    """
    Returns relative response displacement of a linear SDOF system by performing the convolution of a signal and impulse response 
    function, defined in [1]. The function is used in calculation of the extreme response spectrum (ERS) of a random time signal.

    Literature: 
        [1] WILLIAM T. THOMSON, Theory of vibration with applications -> see page 111/512 equation (4.2-5)
    
    :param time_data: signal time data [m/s^2]
    :param dt: time step [s]
    :param f_0: system natural frequency [Hz]
    :param damp: damping ratio [/]

    :return: relative response displacement [m]
    """
    n = len(time_data)
    time = np.arange(n) * dt
    
    omega_0 = 2 * np.pi * f_0
    omega_0d = omega_0 * np.sqrt(1 - damp**2)
    
    impulse_resp_func = -1 / omega_0d * np.exp(-damp * omega_0 * time) * np.sin(omega_0d * time)

    z = signal.convolve(time_data, impulse_resp_func)[:len(time)] * dt
    
    return z


def psd_averaging(self):
    """
    PSD averaging method: Welch's method for calculating PSD of a random signal frm time data.
    """

    if not hasattr(self, 'bins'):
        raise ValueError('Number of bins ``bins`` must be provided for PSD averaging method.')
    
    freq_avg, psd_avg = signal.welch(
        self.time_data, 
        fs=1 / self.dt, 
        nperseg=len(self.time_data) // self.bins, 
        window='boxcar', 
        scaling='density',
        )
    
    self.psd_data = psd_avg
    self.psd_freq = freq_avg

def material_parameters_convert(sigma_f, b, range = False):
    """
    Converts Basquin equation parameters ``sigma_f`` and ``b`` to fatigue life parameters ``C`` and ``k``,
    using a function from FLife package. Basic form of Basquin equation is used here: ``sigma_a = sigma_f* (2*N)**b``. The function converts to parameters from equation ``N * s**k = C``

    :param sigma_f:
        Fatigue strength coefficient [MPa**k].
    :param b:
        Fatigue strength exponent [/]. Represents S-N curve slope.
    :param range:
        False/True sets returned value C with regards to amplitude / range count, respectively.
    
    :return C,k:
        C - S-N curve intercept [MPa**k], k - S-N curve inverse slope [/].

    """

    C,k = basquin_to_sn(sigma_f, b, range=range)
    
    return C, k 
