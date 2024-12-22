..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _CVODE.Alg.Nordsieck:

The Nordsieck History Array
---------------------------

A step from time :math:`t_{n-1}` to :math:`t_n` with step size :math:`h_n =
t_{n} - t_{n-1}` using a method of order :math:`q` begins with the polynomial
interpolant, :math:`\pi_{n-1}`, of degree :math:`q` or less. For Adams methods,
the interpolant satisfies the :math:`q + 1` conditions

.. math::
   :label: cvode_adams_predictor_polynomial

   \pi_{n-1}(t_{n-1}) &= y_{n-1} \\
   \dot{\pi}_{n-1}(t_{n-j}) &= \dot{y}_{n-j}, \quad j = 1,2,\ldots,q.

Similarly for BDF methods, the interpolant satisfies

.. math::
   :label: cvode_bdf_predictor_polynomial

   \pi_{n-1}(t_{n-j}) &= y_{n-j}, \quad j = 1,2,\ldots,q \\
   \dot{\pi}_{n-1}(t_{n-1}) &= \dot{y}_{n-1}.

The solution at :math:`t_n` is obtained by constructing the polynomial
interpolant, :math:`\pi_n`, of degree :math:`q` or less. For Adams methods, the
new interpolant satisfies the :math:`q + 2` conditions

.. math::
   :label: cvode_adams_corrector_polynomial

   \pi_{n}(t_n) &= y_n, \\
   \pi_{n}(t_{n-1}) &= y_{n-1}, \\
   \dot{\pi}_n(t_{n-j}) &= \dot{y}_{n-j} \quad j = 0,1,\ldots,q-1.

Similarly for variable coefficient BDF methods, the new interpolant satisfies

.. math::
   :label: cvode_bdf_corrector_polynomial

   \pi_{n}(t_{n-j}) &= y_{n-j}, \quad j = 0,1,\ldots,q \\
   \dot{\pi}_{n}(t_n) &= \dot{y}_n.

With either method the solution history is represented using the Nordsieck array
given by

.. math::
   :label: cvode_nordsieck_array

   Z_{n-1} = [ y_{n-1},\, h_n \dot{y}_{n-1},\, \frac{1}{2} h_n^2 \ddot{y}_{n-1},\, \ldots,\, \frac{1}{q!} h_n^q y^{(q)}_{n-1}]

where the derivatives, :math:`y^{(j)}_{n-1}` are approximated by
:math:`\pi^{(j)}_{n-1}(t_{n-1})`. With this formulation, computing the time step
is now a problem of constructing :math:`Z_n` from :math:`Z_{n-1}`. This
computation is done by forming the predicted array, :math:`Z_{n(0)}`, with
columns :math:`\frac{1}{j!} h^{(j)}_{n} \pi^{(j)}_{n-1}(t_n)` and computing the
coefficients :math:`l = [l_0,\, l_1,\, \ldots,\, l_q]` such that

.. math::
   :label: cvode_nordsieck_update

   Z_n = Z_{n(0)} + \Delta_n\, l,

where :math:`\Delta_n` is the correction to the predicted solution i.e.,

.. math::
   :label: cvode_correction

   \Delta_n = y_n - y_{n(0)} = \pi_n(t_n) - \pi_{n-1}(t_n)

From the second column of :eq:`cvode_nordsieck_update` we obtain the nonlinear
system for computing :math:`y_n`,

.. math::
   :label: cvode_nonlinear_system

   F(y_n) = y_n - \gamma f(t_n, y_n) - a_n = 0,

where :math:`\gamma = h_n / l_1` and :math:`a_n = y_{n(0)} - \gamma
\dot{y}_{n(0)}`.


Changing the State Size
^^^^^^^^^^^^^^^^^^^^^^^

After completing a step from time :math:`t_{n-1}` to :math:`t_n` using a method
of order :math:`q` and before starting the next step from :math:`t_{n}` to
:math:`t_{n+1}` with a method of order :math:`q'`, the size of the state vector
may change. To "resize" the integrator without restarting from first order, we
require, the user supply, depending on the method, either the recent solution or
right-hand side history for the new state size.

Continuing the integration with the updated state size requires constructing a
new Nordsieck array, :math:`Z_n`, at given the history necessary to build the
appropriate interpolating polynomial, :math:`\pi_n`. To compute :math:`\pi_n`
and its derivatives, we use a Newton interpolating polynomial. Given the
:math:`k+1` data points :math:`x_n,\ldots,x_{n-k}` (solution or right-hand side
values depending on the method) and corresponding times,
:math:`t_n,\ldots,t_{n-k}`, the :math:`k`-th degree polynomial is given by

.. math::
   :label: newton_polynomial

   P(t) = \sum_{j=0}^{k} c_j N_j(t).

The polynomial coefficients, :math:`c_j`, are given by the divided differences,

.. math::
   :label: newton_polynomial_coefficients

   c_j = [x_n,\ldots,x_{n-j}],

where :math:`[x_n,\ldots,x_{n-j}]` is defined recursively as

.. math::
   :label: divided_differences

   [x_i,\ldots,x_{i-j}] = \frac{[x_{i},\ldots , x_{i-j+1}] - [x_{i-1},\ldots , x_{i-j}]}{t_i - t_{i-j}}.

with :math:`[x_i] = x_i`. The basis polynomials, :math:`N_j(t)`, are given by

.. math::
   :label: newton_polynomial_basis

   N_j(t) = \prod_{i=0}^{j-1} (t - t_{n-i}),

for :math:`j > 0` with :math:`N_0(t) = 1`.

Organizing the divided differences in a table illustrates the dependencies in
the recursive computation. For example with :math:`k = 3`, we need to compute
:math:`c_0 = [x_n]`, :math:`c_1 = [x_n,x_{n-1}]`, :math:`c_2 =
[x_n,x_{n-1},x_{n-2}]`, and :math:`c_3 = [x_n,x_{n-1},x_{n-2},x_{n-3}]` which
depend on the difference of the entries to the immediate lower and upper left in
the table below.

.. math::
   :label: divided_differences_table

   \begin{matrix}
   t_n     & x_n     & : & [x_n]     &      &                   &      &                           &      & \\
           &         & : &           & \rhd & [x_n,x_{n-1}]     &      &                           &      & \\
   t_{n-1} & x_{n-1} & : & [x_{n-1}] &      &                   & \rhd & [x_n,x_{n-1},x_{n-2}]     &      & \\
           &         & : &           & \rhd & [x_{n-1},x_{n-2}] &      &                           & \rhd & [x_n,x_{n-1},x_{n-2},x_{n-3}] \\
   t_{n-2} & x_{n-2} & : & [x_{n-2}] &      &                   & \rhd & [x_{n-1},x_{n-2},x_{n-3}] &      & \\
           &         & : &           & \rhd & [x_{n-2},x_{n-3}] &      &                           &      & \\
   t_{n-3} & x_{n-3} & : & [x_{n-3}] &      &                   &      &                           &      &
   \end{matrix}

As such, the coefficients can be computed recursively with the following steps:

1. For :math:`i = 0,\ldots,k`

   :math:`c_i = x_i`

2. For :math:`i = 1,\ldots,k`

   a. For :math:`j = k,\ldots,i`

      :math:`c_j = \frac{c_{j-1} - c_j}{t_{j - i} - t_{j}}`

Rewriting the Newton polynomial using nested multiplications,

.. math::
   :label: newton_nested

   P(t) = c_0 + (t - t_n) \Biggl[ c_1 + (t - t_{n-1}) \biggl[ c_2 + (t - t_{n-2}) \Bigl[ \ldots \bigl[ c_{k-1} + (t - t_{ n - k + 1}) c_k \bigr] \Bigr] \biggr] \Biggr],

leads the following iteration to evaluate :math:`P(t)`:

1. Let :math:`P_k(t) = c_k`

2. For :math:`j = k-1, \ldots, 0`

   :math:`P_{j}(t) = c_j + (t - t_{n-j}) P_{j+1}(t)`

Utilizing this recursive relationship for :math:`P(t)`, we can similarly compute
its derivatives. The evaluation of the :math:`d`-th derivative is given by the
following iteration:

1. Let :math:`P_k(t) = c_k` and :math:`P'_k(t) = P''_k(t) = \ldots = P^{(d)}_k(t) = 0`

2. For :math:`j = k-1, \ldots, 0`

   :math:`P_{j}(t) = c_j + (t - t_{n-j}) P_{j+1}(t)`

   :math:`P'_{j}(t) = P_{j+1}(t) + (t - t_{n-j}) P'_{j+1}(t)`

   :math:`P''_{j}(t) = 2\, P'_{j+1}(t) + (t - t_{n-j}) P''_{j+1}(t)`

   :math:`\vdots`

   :math:`P^{(d)}_{j}(t) = d\, P^{(d-1)}_{j+1}(t) + (t - t_{n-j}) P^{(d)}_{j+1}(t)`

While only :math:`P^{(d)}_k(t) = 0` is needed to start the iteration, note that
:math:`P^{(d)}_{k} = P^{(d)}_{k - 1} = \ldots = P^{(d)}_{k - d + 1} = 0` so some
computations can be skipped when :math:`j > k - d`. With these pieces in place
we can compute the values necessary to build :math:`Z_n`.

For Adams methods we require :math:`y_n` and :math:`q'` right-hand side values,
:math:`\dot{y}_{n-j}` for :math:`j = 0,1,\ldots,q'-1`. The first two columns of
:math:`Z_n` are simply :math:`y_n` and :math:`h_n \dot{y}_n`. To fill the
remaining :math:`q' - 1` columns of :math:`Z_n` we do not need to build
:math:`\pi_n` (as the term interpolating :math:`y_n` vanishes when taking the
first derivative) and instead build :math:`P(t) = \pi'_n` interpolating the
right-hand side values (:math:`x_{n-j} = \dot{y}_{n-j}` above). From
:math:`P(t)` when can then evaluate the higher-order derivatives needed in
:math:`Z_n`.

With BDF methods we need :math:`\dot{y}_n` and :math:`q'` solution values,
:math:`y_{n-j}` for :math:`j = 0,1,\ldots,q'-1`. Again, the first two columns of
:math:`Z_n` are simply :math:`y_n` and :math:`h_n \dot{y}_n`. To fill the
remaining :math:`q' - 1` columns of :math:`Z_n` we need a slight modification of
the Newton polynomial evaluation procedure described above to incorporate the
interpolation condition, :math:`\pi_n(t_n) = \dot{y}_n`. In this case we have
two data points at :math:`t_n` (the state :math:`y_n` and its derivative
:math:`\dot{y}_n`) and the corresponding divided difference is replaced with
:math:`\dot{y}_n`. For example with :math:`q' = 2`, the table of divided
differences is

.. math::
   :label: divided_differences_table_bdf

   \begin{matrix}
   t_n     & y_n     & : & [y_n]     &      &                       &      & \\
           &         &   &           & \rhd & [y_n,y_n] = \dot{y}_n &      & \\
   t_n     & y_n     & : & [y_n]     &      &                       & \rhd & [y_n,y_n,y_{n-1}] \\
           &         &   &           & \rhd & [y_{n},y_{n-1}]       &      & \\
   t_{n-1} & y_{n-1} & : & [y_{n-1}] &      &                       &      &
   \end{matrix}

Other than this adjustment in computing the :math:`c_1` polynomial coefficient,
the iterations for evaluating :math:`P(t)` and :math:`P^{(d)}(t)` to fill the
remaining columns of :math:`Z_n` remain the same.
