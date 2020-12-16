"""
==============================================================
Differential Equation Solver
==============================================================
This module provides a number of analytical methods to solving
ordinary differential equations, like first order, second order,
third order, fourth order, euler's method, runge kutta, and stromer
verlett. It also provides error estimation methods for the
calculations.
-------
"""
__all__ = [
    'solve_ddy',
    'solve_d3y',
    'solve_d4y',
    'euler',
    'runge_kutta',
    'stormer_verlett'
]

import numpy as np


# DEFINING DIFFERENTIAL EQUATION SOLVING TECHNIQUES
# FOR IMPORT TO LIBRARY FRAMEWORK


# 2nd ORDER LINEAR ORDINARY DIFFERENTIAL EQUATIONS
def solve_ddy(function, conditions, delta_t, iterations=1000):
    # FUNCTION SHOULD BE IN THE FORM Y'' = F(Y', Y, X), RETURNING Y DOUBLE PRIME WITH INPUTS Y PRIME, Y, X
    # CONDITIONS SHOULD BE IN THE FORM [Y'0, Y0]
    t = np.linspace(0, delta_t, iterations)

    v = [conditions[0]]
    y = [conditions[1]]
    for i in range(0, len(t) - 1):
        y.append(y[i] + v[i] * (t[i + 1] - t[i]))  # add new to y: existing y plus slope * change in time
        v.append(v[i] + function(v[i], y[i], t[i]) * (
                    t[i + 1] - t[i]))  # add new to slope: existing slope plus 2nd deriv * change in time
    return y


# 3rd ORDER LINEAR ORDINARY DIFFERENTIAL EQUATIONS
def solve_d3y(function, conditions, delta_t, iterations=1000):
    # FUNCTION SHOULD BE IN THE FORM Y''' = F(Y'', Y', Y, X), RETURNING Y TRIPLE PRIME WITH INPUTS Y DOUBLE PRIME, Y PRIME, Y, X
    # CONDITIONS SHOULD BE IN THE FORM [Y''0, Y'0, Y0]
    t = np.linspace(0, delta_t, iterations)
    a = [conditions[0]]
    v = [conditions[1]]
    y = [conditions[2]]
    for i in range(0, len(t) - 1):
        y.append(y[i] + v[i] * (t[i + 1] - t[i]))  # add new to y: existing y plus slope * change in time
        v.append(v[i] + a[i] * (t[i + 1] - t[i]))  # add new to slope: existing slope plus 2nd deriv * change in time
        a.append(a[i] + function(a[i], v[i], y[i]) * (t[i + 1] - t[
            i]))  # add new to acceleration: existing acceleration plus 3rd derivative * change in time
    return y


# 4th ORDER LINEAR ORDINARY DIFFERENTIAL EQUATIONS
def solve_d4y(function, conditions, delta_t, iterations=1000):
    # FUNCTION SHOULD BE IN THE FORM Y'''' = F(Y''', Y'', Y', Y, X), RETURNING 4TH DERIV OF Y WITH INPUTS Y TRIPLE PRIME, Y DOUBLE PRIME, Y PRIME, Y, X
    # CONDITIONS SHOULD BE IN THE FORM [Y'''0, Y''0, Y'0, Y0]
    t = np.linspace(0, delta_t, iterations)
    j = [conditions[0]]
    a = [conditions[1]]
    v = [conditions[2]]
    y = [conditions[3]]
    for i in range(0, len(t) - 1):
        y.append(y[i] + v[i] * (t[i + 1] - t[i]))  # add new to y: existing y plus slope * change in time
        v.append(v[i] + a[i] * (t[i + 1] - t[i]))  # add new to slope: existing slope plus 2nd deriv * change in time
        a.append(a[i] + j[i] * (t[i + 1] - t[
            i]))  # add new to acceleration: existing acceleration plus 3rd derivative * change in time
        j.append(j[i] + function(j[i], a[i], v[i], y[i]) * (t[i + 1] - t[
            i]))  # add new to jerk: existing jerk plus 4th derivative, output of differential function, * change in time
    return y


# EULER's METHOD
def euler(function, y0, delta_t, iterations=1000):
    # FUNCTION SHOULD BE IN THE FORM Y'=F(Y,X), RETURNING Y PRIME WITH INPUTS Y, X
    t = np.linspace(0, delta_t, iterations)
    y = [y0]
    for i in range(0, len(t) - 1):
        y.append(y[i] + function(y[i], t[i]) * (t[i + 1] - t[i]))
    return y


# RUNGE-KUTTA METHOD
def runge_kutta(function, y0, delta_t, iterations=1000):
    # FUNCTION SHOULD BE IN THE FORM Y'=F(Y,X), RETURNING Y PRIME WITH INPUTS Y, X
    t = np.linspace(0, delta_t, iterations)
    y = [y0]
    for i in range(0, len(t) - 1):
        h = t[i + 1] - t[i]
        k1 = function(y[i], t[i])
        k2 = function(y[i] + k1 * h / 2, t[i] + h / 2)
        k3 = function(y[i] + k2 * h / 2, t[i] + h / 2)
        k4 = function(y[i] + k3 * h, t[i] + h)

        y.append(y[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4))
    return y


# graph test runge-kutta


# STORMER VERLETT METHOD
def stormer_verlett(function, conditions, delta_t, iterations=1000):
    # FUNCTION SHOULD BE IN THE FORM Y'' = F(Y), RETURNING Y DOUBLE PRIME WITH INPUT Y
    # CONDITIONS SHOULD BE IN THE FORM [Y'0, Y0]
    t = np.linspace(0, delta_t, iterations)

    v = [conditions[0]]
    y = [conditions[1]]
    y.append(y[0] + v[0] * (t[1] - t[0]) + 0.5 * function(y[0]))  # FIRST iterations
    for i in range(1, len(t) - 1):
        y.append(2 * y[i] - y[i - 1] + function(y[i]) * (t[1] - t[0]) ** 2)
    return y

