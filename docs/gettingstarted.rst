.. _gettingstarted:

=================
 Getting started
=================

The QEH python package implements a simple command line
interface (CLI) to make QEH calculations without writing python
scripts. This is very useful when investigating standard properties of
heterostructures. However, QEH calculations can also be performed in
Python scripts which is useful if the CLI is not sufficient. This
tutorial first introduces the QEH CLI and ends with a description of
writing Python scripts that call the QEH code.


Command line interface
======================
To see the help for the CLI do

.. command-output:: qeh -h
   :ellipsis: 10

There are many options for the CLI but most important is the
``layers`` argument which describes the layers that make up the
heterostructure. For example to build a heterostructure consisting of
doped graphene on 3 layers of boron nitride do::

  qeh graphene+doping=0.5 3BN

For example the plasmons of a doped graphene on three layers of boron
nitride can be calculated and saved as::

  qeh graphene+doping=0.5 3BN --plasmons --save-plots "_g-3BN"

which creates a figure ``loss_g-3BN.png`` containing the loss spectrum
and a figure ``plasmon_modes_g-3BN.png`` which traces the plasmon
modes:

.. image:: loss_g-3BN.png
   :width: 49%
	
.. image:: plasmon_modes_g-3BN.png
   :width: 49%

Here we introduced a couple of key concepts. You can modify layers
with the ``+`` syntax. For example, the example above used the layer
modifier ``doping=0.5`` which doped graphene. See the doping_ section
for details about this modifier. Because boron nitride is polar there
will be an additional contribution to the dielectric properties from
its phonons. These have been calculated as well and can be included
with the phonons_ modifier

Available monolayers
====================
Here you find a list of the monolayers that are available for the QEH
model: XXX

Doping
======
.. _doping:
The doping modifier can be used to dope a layer. The doping modifier
take additional arguments and an example could be::


  

Phonons
=======
.. _phonons:



Doped semiconductors
====================
The QEH package contain effective masses for all TMDCs which can be
used in connection with the qeh model to calculate the properties of
doped semi-conductors::

  qeh H-MoS2+doping=0.1,T=25e-3,eta=1e-3 --plasmons --save-fig '_doped_MoS2'

Here we dope MoS2 in the H phase with at a finite temperature
``T=25e-3`` of 25 meV and a relaxation rate of ``eta=1e-3`` 3 meV.
