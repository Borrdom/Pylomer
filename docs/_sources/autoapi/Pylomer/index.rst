:py:mod:`Pylomer`
=================

.. py:module:: Pylomer


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   FloryHuggins/index.rst
   GordonTaylor/index.rst
   crank_and_hopfenberg/index.rst


Package Contents
----------------


Functions
~~~~~~~~~

.. autoapisummary::

   Pylomer.crank
   Pylomer.BH_Model
   Pylomer.lnai
   Pylomer.TgGT



.. py:function:: crank(t, D, L, w0=0, winfty=1)

   Model solvent diffusion in a 1D slab according to cranks equation p.48 (equation 4.18)
   Crank, The mathematics of diffusion, 2nd Ed., Clarendon, Oxford, United Kingdom, 1976, 85 (ISBN 0 19 853344 6)

   :param t: time /s
   :type t: array_like
   :param D: Solvent diffusion coefficient  / m^2/s
   :type D: float
   :param L: Thickness /m
   :type L: float
   :param w0: solvent mass fraction at t=0
   :type w0: float
   :param winfty: solvent mass fraction in equilibrium
   :type winfty: float

   :returns: solvent mass fraction at t /
   :rtype: ndarray


.. py:function:: BH_Model(t, D, L, kr, difffrac, w0=0, winfty=1)

   Model solvent diffusion in a 1D slab with Hophenberg modification for relaxation https://doi.org/10.1016/0032-3861(78)90269-0

   :param t: time /s
   :type t: array_like
   :param D: Solvent diffusion coefficient  / m^2/s
   :type D: float
   :param L: Thickness /m
   :type L: float
   :param kr: time constant for relaxation
   :type kr: float
   :param difffrac: specifies the fractional solvent amount that is sorbed soley through diffusion
   :type difffrac: float
   :param w0: solvent mass fraction at t=0
   :type w0: float
   :param winfty: solvent mass fraction in equilibrium
   :type winfty: float

   :returns: solvent mass fraction at t /
   :rtype: ndarray


.. py:function:: lnai(wi, rho0i, Mi, chi=None)

   Compute the log activity of a mixture

   :param wi: weight fractions number of components
   :type wi: array_like
   :param Mi: Molar mass / g/mol
   :type Mi: array_like
   :param chi: Flory Huggins Chi Prameter        /-
   :type chi: optional,array_like

   :returns: logarithmic activity of component i  /K
   :rtype: ndarray


.. py:function:: TgGT(wi, Tg0i, q=0, Ki=None, rho0i=None)

   Compute the glass transition temperature of a mixture

   :param wi: 2D Array of weight fractions [ number of components,number of Points]
   :type wi: array_like
   :param Tg0i: pure component glass transition temperature /K
   :type Tg0i: array_like
   :param q: Kwei parameter /-
   :type q: array_like
   :param rho0i: pure component densities /kg/m^3
   :type rho0i: optional,array_like
   :param Ki: Gordon-Taylor parameters         /-
   :type Ki: optional,array_like

   :returns: glass transition temperature of a mixture  /K
   :rtype: ndarray


