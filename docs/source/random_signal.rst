Random signal
================

Here is an example of determining the ERS and the FDS of a random signal, defined in the time and in the frequency domain.

Import the required packages
----------------------------


.. code-block:: python

    import numpy as np
    import pyFDS
    import matplotlib.pyplot as plt
    import pyExSi as es

SpecificationDevelopment object
-------------------------------

Specification development object contains data, required for the calculation of extreme response spectrum (ERS) and fatigue damage spectrum (FDS). It enables calculation for random signals, that are defined with the PSD or with the time-history.
For the time-history, two methods are available. ERS and FDS can be determined directly from time history using convolution, or by first converting the history into PSD and calculating spectra from the PSD.

This example contains an ERS and FDS calculation with all three avaiable methods:

* from PSD

* from time-history using convolution (directly from time history)

* from time-history using PSD averaging (conversion to PSD, then ERS and FDS from PSD)

All of the ERS and FDS are plotted so the methods can be compared.

Random signal (flat-shaped PSD)
--------------------------------

Generate data
~~~~~~~~~~~~~

Random signal is generated using PyExSi.

.. code-block:: python

    fs = 5000 # sampling frequency [Hz]
    time= 500

    N = int(time*fs) # number of data points of time signal
    t = np.arange(0,N)/fs # time vector

    # define frequency vector and one-sided flat-shaped PSD
    freq_flat = np.arange(0, fs/2, 1/time) # frequency vector
    freq_lower = 200 # PSD lower frequency limit  [Hz]
    freq_upper = 1000 # PSD upper frequency limit [Hz]
    PSD_flat = es.get_psd(freq_flat, freq_lower, freq_upper,variance=800) # one-sided flat-shaped PSD

    #get gaussian stationary signal
    gausian_signal = es.random_gaussian(N, PSD_flat, fs)

Plot the generated signal:

.. code-block:: python

    plt.plot(freq_flat, PSD_flat)
    plt.xlabel('f [Hz]')
    plt.ylabel('(m/s²)²/Hz')
    plt.show()

    plt.plot(t,gausian_signal)
    plt.xlabel('t [s]')
    plt.ylabel('m/s²')
    plt.show()

Instantiate the SpecificationDevelopment object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Object is instantiated with inputs:

* ``freq_data``: tuple containing (``f0_start``, ``f0_stop``, ``f0_step``) [Hz] or a frequency vector (array), defining the range where the ERS and FDS will be calculated.

*  damping ratio ``damp`` or damping Q-factor ``Q``.

.. code-block:: python
    
    sd_flat_psd = pyFDS.SpecificationDevelopment(freq_data=(100,1100,20),damp=0.05)
    sd_flat_time = pyFDS.SpecificationDevelopment(freq_data=(100,1100,20),damp=0.05)
    sd_flat_time_2 = pyFDS.SpecificationDevelopment(freq_data=(100,1100,20),damp=0.05)

Set the random load
~~~~~~~~~~~~~~~~~~~

Random load is defined with the ``set_random_load`` method. Time history or PSD must given as input. Class method automatically determines, whether the input is time history or PSD, based on the type of input:

* PSD: input is tuple containing (psd data (array), frequency vector (array)).

* Time history: input is tuple containing (time history data (array), dt (scalar)).

If time history is given as input, method of spectra calculation must also be defined. Available methods are:

* ``convolution`` (directly from time history)

* ``psd_averaging`` (conversion to PSD, then ERS and FDS from PSD)

.. code-block:: python

    sd_flat_psd.set_random_load((PSD_flat[::100],freq_flat[::100]),unit='ms2',T=500) #input is tuple (psd array, freq array)
    sd_flat_time.set_random_load((gausian_signal,1/fs), unit='ms2',method='convolution') #input is tuple (psd data, frequency vector)
    sd_flat_time_2.set_random_load((gausian_signal,1/fs), unit='ms2',method='psd_averaging',bins=500) #input is tuple (psd data, frequency vector)

Get the ERS and FDS
~~~~~~~~~~~~~~~~~~~~

ERS and FDS are calculated with the ``get_ers`` and ``get_fds`` methods. For the FDS calculation, the additional material fatigue parameters ``b``, ``C`` and ``K`` must be provided.

.. code-block:: python
    
    sd_flat_psd.get_ers()
    sd_flat_time.get_ers()
    sd_flat_time_2.get_ers()

    b=10
    C=1e80
    K=6.3*1e10

    sd_flat_psd.get_fds(b=b,C=C,K=K)
    sd_flat_time.get_fds(b=b,C=C,K=K)
    sd_flat_time_2.get_fds(b=b,C=C,K=K)

Plot the results
~~~~~~~~~~~~~~~~

ERS and FDS are plotted for all three methods.

.. code-block:: python

    plt.plot(sd_flat_psd.f0_range,sd_flat_psd.ers,label='PSD')
    plt.plot(sd_flat_time.f0_range,sd_flat_time.ers,label='Time history (convolution)')
    plt.plot(sd_flat_time_2.f0_range,sd_flat_time_2.ers,label='Time history (psd averaging)')
    plt.title('ERS')
    plt.legend()
    plt.grid()
    plt.ylabel('[m/s²]')
    plt.xlabel('f [Hz]')
    plt.show()

    plt.semilogy(sd_flat_psd.f0_range,sd_flat_psd.fds,label='PSD')
    plt.semilogy(sd_flat_time.f0_range,sd_flat_time.fds,label='Time history (convolution)')
    plt.semilogy(sd_flat_time_2.f0_range,sd_flat_time_2.fds,label='Time history (psd averaging)')
    plt.title('FDS')
    plt.ylabel('Damage')
    plt.xlabel('f [Hz]')
    plt.grid()

