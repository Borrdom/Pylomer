:py:mod:`Pylomer.crank_and_hopfenberg`
======================================

.. py:module:: Pylomer.crank_and_hopfenberg


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   Pylomer.crank_and_hopfenberg.crank
   Pylomer.crank_and_hopfenberg.BH_Model



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


