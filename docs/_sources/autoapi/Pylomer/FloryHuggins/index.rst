:py:mod:`Pylomer.FloryHuggins`
==============================

.. py:module:: Pylomer.FloryHuggins


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   Pylomer.FloryHuggins.lngi



.. py:function:: lngi(wi, Mi, rho0i, chi=None)

   Compute the glass transition temperature of a mixture

   :param wi: 2D Array of weight fractions [number of components,number of Points]
   :type wi: array_like
   :param Mi: pure component glass transition temperature /K
   :type Mi: array_like
   :param rho0i: pure component densities /kg/m^3
   :type rho0i: optional,array_like
   :param Ki: Gordon-Taylor parameters         /-
   :type Ki: optional,array_like

   :returns: glass transition temperature of a mixture  /K
   :rtype: ndarray


