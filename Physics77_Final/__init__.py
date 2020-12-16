"""
Physics77_Final
=====
Provides
  1. Differential Equation Library(Ordered Equations, Euler's Method, Runge Kutta, Stormer Verlet)
  2. Multibody systems and graphing
  3. Lagrange Point Solver, Error Estimation
How to use the documentation
----------------------------
Documentation is available in docstrings provided.
We recommend exploring the docstrings using
`IPython <https://ipython.org>`_, an advanced Python shell with
TAB-completion and introspection capabilities.  See below for further
instructions.
The docstring examples assume that `Physics77_Final` has been imported as `pf`::
  >>> import Physics77_Final as pf
Code snippets are indicated by three greater-than signs::
  >>> x = 42
  >>> x = x + 1
Use the built-in ``help`` function to view a function's docstring::
  >>> help(pf.coupled_equation_solver_1D)
  ... # doctest: +SKIP
Available subpackages/modules
---------------------
diffeq
    Analytical solvers to commonly found differential equations including
    first and second order equations, Euler's method, Runge Kutta Method,
    Stromer Verlet Method. Error estimate functions for these methods are
    also provided.
multibody
    Graphs and animates out solutions to multibody systems
---------------------
Viewing documentation using IPython
"""