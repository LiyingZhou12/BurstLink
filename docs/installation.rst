Installation
============

Please follow the guide below to install BurstLink and its dependent
software.

System requirements
~~~~~~~~~~~~~~~~~~~

BurstLink was developed and tested on Linux (U), macOS and Windows
operating systems. Linux and macOS is recommended.

Python requirements
~~~~~~~~~~~~~~~~~~~

BurstLink was developed and tested using python 3.8. We recommend using
mamba or conda to manage Python environments.

Installation using ``pip``
~~~~~~~~~~~~~~~~~~~~~~~~~~

We suggest setting up BurstLink in a separate ``mamba`` or ``conda``
environment to prevent conflicts with other software dependencies.
Create a new Python environment specifically for BurstLink and install
the required libraries within it.

.. code:: bash

   mamba create -n burstlink_env python=3.8 r-base=4.3.2
   mamba activate burstlink_env
   pip install burstlink

if you use ``conda``, ``r-base=4.3.2`` may not included in the channels.
Instead, you can ``r-base=4.3.1`` in ``conda``.

Now you can check if burstlink can sucessfully import in Python:

.. code:: python

   import burstlink as bl

Troubleshooting
~~~~~~~~~~~~~~~

**Continuous updating…**

Feel free to contact us if you encounter any other problems with the
installation.
