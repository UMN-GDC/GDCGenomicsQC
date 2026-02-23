.. highlight:: shell

============
Installation
============


Stable release
--------------

When the software is complete, run this command in your terminal:

.. code-block:: console

    $ pip install GDC


From Source (Development)
-------------------------

If you wish to contribute or use the latest developmental features, install directly from the source code.

1. **Clone the repository:**

   .. code-block:: console

       $ git clone https://github.com/and02709/GDC.git
       $ cd GDC

2. **Set up a Virtual Environment (Recommended):**

   .. code-block:: console

    $ source /home/gdc/public/envs/load_miniconda3.sh
    $ module load plink
    $ module load perl


External Dependencies
---------------------

GDC acts as a high-level interface for genetic data processing. **You must have the following binaries installed and available in your system PATH.**

PLINK (1.9 and 2.0)
~~~~~~~~~~~~~~~~~~~

GDC relies heavily on PLINK for rapid data cleaning.

* **Download:** Get the latest stable binaries from `cog-genomics.org <https://www.cog-genomics.org/plink/>`_.
* **Setup:** Ensure the executables are named ``plink`` (for 1.9) and ``plink2`` (for 2.0).
* **Verify:** Run the following in your terminal to ensure they are accessible:

  .. code-block:: console

      $ plink --version
      $ plink2 --version

.. note::
   If you are working on a High-Performance Computing (HPC) cluster using SLURM, you may need to load these via modules (e.g., ``module load plink/1.9``).


Configuration & Troubleshooting
-------------------------------

If GDC cannot find your PLINK binaries, you can verify your system PATH:

.. code-block:: console

    # Linux/macOS
    $ echo $PATH

    # Windows (PowerShell)
    $ $env:Path

Ensure the directory containing your binaries is listed. If you encounter permission errors when running from source, ensure the scripts have execution privileges:

.. code-block:: console

    $ chmod +x GDC/*.py
