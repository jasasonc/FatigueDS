Sine and sine sweep signal
===========================

Here are the examplea of determining the ERS and FDS of a sine and sine-sweep signal:

Sine signal
------------

Import the required packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. code-block:: python

    import numpy as np
    import pyFDS
    import matplotlib.pyplot as plt


Instantiate SpecificationDevelopment object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SpecificationDevelopment object contains data, required for the calculation of extreme response spectrum (ERS) and fatigue damage spectrum (FDS).

.. code-block:: python

    sd_sine = pyFDS.SpecificationDevelopment(freq_data=(0,2000,5), damp=0.1)

Set the sine load
~~~~~~~~~~~~~~~~~~

Set the sine load with the following parameters:

* ``sine_freq`` - frequency of the sine signal

* ``amp`` - amplitude of the sine signal

* ``t_total`` - total time of the signal

* ``exc_type`` - excitation type

.. code-block:: python

    sd_sine.set_sine_load(sine_freq=500,amp=10,t_total=3600, exc_type='acc')

Calculate the ERS and FDS
~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate the ERS and FDS. For the FDS calculation, the additional mateirla fatigue parameters ``b``, ``C`` and ``K`` must be provided.

.. code-block:: python

    sd_sine.get_ers()
    sd_sine.get_fds(b=10,C=1e80,K=6.3*1e10)

Plot the results
~~~~~~~~~~~~~~~~~

Plot the ERS and FDS results.

.. code-block:: python

    fig, axs = plt.subplots(2,1, figsize=(6, 8))
    axs[0].plot(sd_sine.f0_range, sd_sine.ers)
    axs[0].set_ylabel('ERS [m/s^2]')

    axs[1].loglog(sd_sine.f0_range, sd_sine.fds)
    axs[1].set_ylabel('FDS [damage]')


Sine-sweep signal
------------------

Instantiate SpecificationDevelopment object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SpecificationDevelopment object contains data, required for the calculation of extreme response spectrum (ERS) and fatigue damage spectrum (FDS).

.. code-block:: python

    sd_sine_sweep = pyFDS.SpecificationDevelopment(freq_data=(0,2000,5), damp=0.1)


Set the sine-sweep load
~~~~~~~~~~~~~~~~~~~~~~~~

Set the sine-sweep load with the following parameters:

* ``const_amp`` - amplitude of the sine-sweep signal

* ``const_f_range`` - frequency range of the sine-sweep signal

* ``exc_type`` - excitation type

* ``dt`` - time step

* ``sweep_type`` - sweep type

* ``sweep_rate`` - sweep rate

.. code-block:: python

    sd_sine_sweep.set_sine_sweep_load(const_amp=[5,10,20], const_f_range=[20,100,500,1000],exc_type='acc', dt=1, sweep_type='log', sweep_rate=1)


Calculate the ERS and FDS
~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate the ERS and FDS. For the FDS calculation, the additional mateirla fatigue parameters ``b``, ``C`` and ``K`` must be provided.

.. code-block:: python

    sd_sine_sweep.get_ers()
    sd_sine_sweep.get_fds(b=10,C=1e80,K=6.3*1e10)

Plot the results
~~~~~~~~~~~~~~~~~

Plot the ERS and FDS results.

.. code-block:: python

    fig, axs = plt.subplots(2,1, figsize=(6, 8))
    axs[0].plot(sd_sine_sweep.f0_range, sd_sine_sweep.ers)
    axs[0].set_ylabel('ERS [m/s^2]')

    axs[1].loglog(sd_sine_sweep.f0_range, sd_sine_sweep.fds)
    axs[1].set_ylabel('FDS [Damage]')

    plt.show()