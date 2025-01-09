pyFDS
-----------------------

Calculating Extreme Response Spectrum (ERS) and Fatigue Damage Spectrum (FDS) of signals. 
Calculation supported for sine, sine-sweep and random signals (defined with PSD or time history).
Theory based on [1].


Installation
------------------

Use `pip` to install it by:

.. code-block:: console

    $ pip install pyFDS

Usage
------------------
Some short examples of how to use the package are given below for different types of signals.

Random signals (PSD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is an example of determining the ERS and FDS of a random signal, defined in the frequency domain (PSD):

.. code-block:: python

    import numpy as np
    import pyFDS
    import pyExSi as es
    import matplotlib.pyplot as plt

    # generate random signal

    fs = 5000 # sampling frequency [Hz]
    time= 1 # time duration [s]

    # define frequency vector and one-sided flat-shaped PSD
    freq_flat = np.arange(0, fs/2, 1/time) # frequency vector
    freq_lower = 200 # PSD lower frequency limit  [Hz]
    freq_upper = 1000 # PSD upper frequency limit [Hz]
    PSD_flat = es.get_psd(freq_flat, freq_lower, freq_upper,variance=800) # one-sided flat-shaped PSD

    # instantiate the SpecificationDevelopment class
    sd_flat_psd = pyFDS.SpecificationDevelopment(freq_data=(100,1100,20),damp=0.05)

    # set the random load
    sd_flat_psd.set_random_load((PSD_flat,freq_flat),unit='ms2',T=3600) # input is PSD and frequency vector

    # calculate the ERS and FDS
    sd_flat_psd.get_ers()
    sd_flat_psd.get_fds(b=10,C=1e80,K=6.3*1e10)
    
    #plot the results
    plt.plot(sd_flat_psd.f0_range,sd_flat_psd.ers)
    plt.title('ERS')
    plt.grid()
    plt.ylabel('[m/s²]')
    plt.xlabel('f [Hz]')
    plt.show()

    plt.semilogy(sd_flat_psd.f0_range,sd_flat_psd.fds,label='PSD')
    plt.title('FDS')
    plt.ylabel('Damage')
    plt.xlabel('f [Hz]')
    plt.grid()

Random signals (time history)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is an example of determining the ERS and FDS of a random signal, defined in the time domain. For time domain, two methods are available:
    - Convolution (directly from time history, using rainflow counting)
    - PSD averaging (converting time history to PSD and then to ERS and FDS)

.. code-block:: python

    import numpy as np
    import pyFDS
    import matplotlib.pyplot as plt

    # import random signal
    _time_data = np.load('test_data/test_time_history.npy', allow_pickle=True)
    time_history_data = _time_data[:,1]
    t = _time_data[:,0] 
    dt = t[2] - t[1]

    #instantiate the SpecificationDevelopment classes
    sd_1 = pyFDS.SpecificationDevelopment(freq_data=(20,200,5)) #convolution
    sd_2 = pyFDS.SpecificationDevelopment(freq_data=(20,200,5)) #psd averaging

    # set the random loads (input is time history data and time step)
    sd_1.set_random_load((time_history_data,dt), unit='g', method='convolution')
    sd_2.set_random_load((time_history_data,dt), unit='g',method='psd_averaging',bins=10)

    # calculate the ERS and FDS
    sd_1.get_ers()
    sd_1.get_fds(b=10,C=1e80,K=6.3*1e10)

    sd_2.get_ers()
    sd_2.get_fds(b=10,C=1e80,K=6.3*1e10)

    # plot the results
    plt.plot(sd_1.f0_range,sd_1.ers,label='Time history (convolution)')
    plt.plot(sd_2.f0_range,sd_2.ers,label='Time history (PSD averaging)')
    plt.title('ERS')
    plt.legend()
    plt.grid()
    plt.ylabel('[g]')
    plt.xlabel('f [Hz]')
    plt.show()

    plt.loglog(sd_1.f0_range,sd_1.fds,label='Time history (convolution)')
    plt.loglog(sd_2.f0_range,sd_2.fds,label='Time history (PSD averaging)')
    plt.title('FDS')
    plt.ylabel('Damage')
    plt.xlabel('f [Hz]')
    plt.grid()
    plt.legend()

Sine and sine-sweep signals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is an example of determining the ERS and FDS of a sine and sine-sweep signal:

.. code-block:: python

    import numpy as np
    import pyFDS
    import matplotlib.pyplot as plt

    #instantiate classes
    sd_sine = pyFDS.SpecificationDevelopment(freq_data=(0,2000,5), damp=0.1) #sine
    sd_sine_sweep = pyFDS.SpecificationDevelopment(freq_data=(0,2000,5), damp=0.1) #sine sweep

    # set the sine and sine-sweep loads
    sd_sine.set_sine_load(sine_freq=500,amp=10,t_total=3600) # t_total is only needed for fds calculation
    sd_sine_sweep.set_sine_sweep_load(const_amp=[5,10,20], const_f_range=[20,100,500,1000],exc_type='acc', sweep_type='log', sweep_rate=1)

    # calculate the ERS and FDS
    sd_sine.get_ers()
    sd_sine_sweep.get_ers()

    sd_sine.get_fds(b=10,C=1e80,K=6.3*1e10)
    sd_sine_sweep.get_fds(b=10,C=1e80,K=6.3*1e10)

    # plot the results
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    axs[0, 0].plot(sd_sine.f0_range, sd_sine.ers)
    axs[0, 0].set_title('Sine')
    axs[0, 0].set_ylabel('ERS [m/s^2]')

    axs[1, 0].loglog(sd_sine.f0_range, sd_sine.fds)
    axs[1, 0].set_ylabel('FDS')
    axs[1, 0].set_xlabel('Frequency Range')

    axs[0, 1].plot(sd_sine_sweep.f0_range, sd_sine_sweep.ers)
    axs[0, 1].set_title('Sine sweep')
    axs[1, 1].loglog(sd_sine_sweep.f0_range, sd_sine_sweep.fds)
    axs[1, 1].set_xlabel('Frequency Range')

    plt.tight_layout()
    plt.show()


References:
    1. C. Lalanne, Mechanical Vibration and Shock: Specification development,
    London, England: ISTE Ltd and John Wiley & Sons, 2009
