Getting started
===============

To install pyFDS, use ``pip``:

.. code-block:: console

    $ pip install pyFDS


Import the package
-------------------

.. code-block:: python
    
    import pyFDS


Create a SpecificationDevelopment object
-----------------------------------------

To calculate the Extreme Response Spectrum (ERS) and Fatigue Damage Spectrum (FDS), SpecificationDevelopment object must be created. The object is created by providing the frequency range and damping ratio.
Frequency range ``freq_data`` is defined by a tuple (``f0_start``, ``f0_stop``, ``f0_step``) in Hz and sets the frequency points where the ERS and FDS will be calculated. Alternatevely, frequency vector (array) can be passed as input. 
Damping ratio ``damp`` is a float value between 0 and 1.

.. code-block:: python

    sd = pyFDS.SpecificationDevelopment(freq_data=(f0_start,f0_stop,f0_step),damp)

    # or

    sd = pyFDS.SpecificationDevelopment(freq_vector,damp)


Setting the load signal
------------------------

Load signal is defined with different methods, depending on the type of signal.

Sine signal
~~~~~~~~~~~~

For the sine signal; frequency, amplitude, total time of the signal and excitation type must be provided.

.. code-block:: python
    
    sd.set_sine_load(sine_freq,amp,t_total, exc_type)


Sine-sweep signal
~~~~~~~~~~~~~~~~~~

For the sine-sweep signal; amplitude, frequency range, excitation type, time step, sweep type and sweep rate must be provided.

.. code-block:: python

    sd.set_sine_sweep_load(const_amp, const_f_range, exc_type, dt, sweep_type, sweep_rate)


Random signal (PSD)
~~~~~~~~~~~~~~~~~~~~

For the random signal defined by the Power Spectral Density (PSD); PSD, frequency range, unit and time duration must be provided.

.. code-block:: python

    sd.set_random_load((PSD,freq), unit, T)


Random signal (time history)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the random signal defined by the time history; time history, time step, unit and method of calculation must be provided.
Available methods are:

* ``convolution`` (Directly from time history)

* ``psd_averaging`` (Conversion to PSD, then to ERS and FDS from PSD)


.. code-block:: python

    sd.set_random_load((time_history,dt), unit, method)


Calculating the ERS and FDS
----------------------------

After the load signal is set, the ERS and FDS can be calculated.

ERS is calculated by:

.. code-block:: python

    sd.get_ers()


FDS calculation requires additional material fatigue parameters: ``b``, ``C`` and ``K``. It is calculated by:

.. code-block:: python

    sd.get_fds(b,C,K)