Sine and sine-sweep signal
===========================

Here are examples of how to determine the ERS and FDS of sine and sine-sweep signals:

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

    sd_1 = pyFDS.SpecificationDevelopment(freq_data=(0, 2000, 5), damp=0.1)

Set the sine load
~~~~~~~~~~~~~~~~~~

Set the sine load with the following parameters:

* ``sine_freq`` - frequency of the sine signal

* ``amp`` - amplitude of the sine signal

* ``t_total`` - total time of the signal

* ``exc_type`` - excitation type. Possible values are ``acc``, ``vel`` and ``disp`` 

.. code-block:: python

    sd_1.set_sine_load(sine_freq=500, amp=10, t_total=3600, exc_type='acc')

Calculate the ERS and FDS
~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate the ERS and FDS. For the FDS calculation, the additional mateiral fatigue parameters ``b``, ``C`` and ``K`` must be provided.

.. code-block:: python

    sd_1.get_ers()
    sd_1.get_fds(b=10, C=1e80, K=6.3 * 1e10)

Plot the results
~~~~~~~~~~~~~~~~~

Plot the ERS and FDS results.

.. code-block:: python

    sd_1.plot_ers(label='sine')
    sd_1.plot_fds(label='sine')

Or access the results directly:

.. code-block:: python

    sd_1.ers
    sd_1.fds
    sd_1.f0_range  # frequency vector


Sine-sweep signal
------------------

Instantiate SpecificationDevelopment object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SpecificationDevelopment object contains data, required for the calculation of extreme response spectrum (ERS) and fatigue damage spectrum (FDS).

.. code-block:: python

    sd_2 = pyFDS.SpecificationDevelopment(freq_data=(0, 2000, 5), damp=0.1)


Set the sine-sweep load
~~~~~~~~~~~~~~~~~~~~~~~~

Set the sine-sweep load with the following parameters:

* ``const_amp`` - amplitude of the sine-sweep signal

* ``const_f_range`` - frequency range of the sine-sweep signal

* ``exc_type`` - excitation type. Possible values are ``acc``, ``vel`` and ``disp`` 

* ``dt`` - time step

* ``sweep_type`` - sweep type. Possible values are ``log`` and ``lin``

* ``sweep_rate`` - sweep rate

.. code-block:: python

    sd_2.set_sine_sweep_load(const_amp=[5, 10, 20], const_f_range=[20, 100, 500, 1000], exc_type='acc', dt=1, sweep_type='log', sweep_rate=1)


Calculate the ERS and FDS
~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate the ERS and FDS. For the FDS calculation, the additional mateirla fatigue parameters ``b``, ``C`` and ``K`` must be provided.

.. code-block:: python

    sd_2.get_ers()
    sd_2.get_fds(b=10, C=1e80, K=6.3 * 1e10)

Plot the results
~~~~~~~~~~~~~~~~~

Plot the ERS and FDS results.

.. code-block:: python

    sd_2.plot_ers(label='sine sweep')
    sd_2.plot_fds(label='sine sweep')

    plt.show()

Or access the results directly:

.. code-block:: python

    sd_2.ers
    sd_2.fds
    sd_2.f0_range  # frequency vector