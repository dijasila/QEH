.. Quantum Electrostatic Heterostructure Model documentation master file, created by
   sphinx-quickstart on Thu Sep 12 08:48:51 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Quantum Electrostatic Heterostructure Model's documentation!
=======================================================================
The quantum electrostatic heterostructure (QEH) model is a fast
and accurate model for simulating optical properties of stacks
of 2D materials (heterostructures). This python package implements
an easy command line interface to quickly calculate hetorostructure
properties.

Minimal example of calculating the plasmons of a graphene boron nitride
heterostructure::

   $ qeh graphene+doping=0.5 3BN graphene+doping=0.5 --plasmons --plot

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install.rst
   tutorials.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
