:py:mod:`Pylomer.FloryHuggins`
==============================

.. py:module:: Pylomer.FloryHuggins


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   Pylomer.FloryHuggins.lnai



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


