:py:mod:`Pylomer.GordonTaylor`
==============================

.. py:module:: Pylomer.GordonTaylor


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   Pylomer.GordonTaylor.TgGT



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


