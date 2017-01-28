.. _home:

kdotp-symmetry
==============

This is a tool to calculate the general form of a :math:`\mathbf{k}\cdot\mathbf{p}` Hamiltonian under a given symmetry constraint.

.. contents ::
    :local:

Usage
-----

Installation
~~~~~~~~~~~~

You can install this tool with with pip:

.. code ::

    pip install kdotp-symmetry
    
Example: TaAs\ :sub:`2`
~~~~~~~~~~~~~~~~~~~~~~~

Its usage is best explained with an example -- I'll pick a four-band Hamiltonian for TaAs\ :sub:`2`. If you're interested in the physics of this particular example, it comes from a `paper <https://arxiv.org/abs/1611.01858>`_ of ours where we also used this :math:`\mathbf{k}\cdot\mathbf{p}` model. As for now, all we will need to know about the material is its symmetry, and the relevant representations for the given bands.

The symmetry group of the material is :math:`C2 / m` (space group 12), which means it has rotation :math:`C_{2y}`, parity :math:`P`, mirror :math:`M_y` and time-reversal symmetry :math:`\mathcal{T}`. For the analysis of the Hamiltonian we only need a generating set of the group, so we can pick :math:`C_{2y}`, :math:`P` and :math:`\mathcal{T}`. In the particular basis we chose for this analysis, the :math:`\mathbf{k}`-space matrix for these symmetries are as follows:

.. math :: 
    
    C_{2y} =& ~\begin{pmatrix} 0&1&0 \\ 1&0&0 \\ 0&0&-1 \end{pmatrix}\\
    P =& ~-\mathbb{1}_{3\times 3}\\
    \mathcal{T} =& ~-\mathbb{1}_{3\times 3}
    
The corresponding representations are 

.. math ::
    
    C_{2y} =& ~\begin{pmatrix} i&0&0&0 \\ 0&-i&0&0 \\ 0&0&i&0 \\ 0&0&0&-i \end{pmatrix} \\
    P =& ~\begin{pmatrix} 1&0&0&0 \\ 0&1&0&0 \\ 0&0&-1&0 \\ 0&0&0&-1 \end{pmatrix} \\
    \mathcal{T} =& ~\begin{pmatrix} 0&-1&0&0 \\ 1&0&0&0 \\ 0&0&0&-1 \\ 0&0&1&0 \end{pmatrix} ~\hat{K},

where :math:`\hat{K}` is complex conjugation.

Creating the :class:`.SymmetryOperation`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to run the code, we must first specify the symmetries as described above. To do this, we can use the :class:`.SymmetryOperation` class, which is a :py:func:`namedtuple <collections.namedtuple>` with two attributes ``kmatrix`` and ``repr``.

The first attribute, ``kmatrix``, is just the :math:`\mathbf{k}`-space matrix for the symmetry. 

The second, ``repr``, describes the symmetry representation which can either be a unitary matrix :math:`U`, or a unitary matrix and complex conjugation :math:`U \hat{K}`. As such, the ``repr`` is another :py:func:`namedtuple <collections.namedtuple>` called :class:`.Representation` with two attributes ``matrix`` and ``complex_conjugate``. The ``matrix`` is the unitary matrix :math:`U`, and ``complex_conjugate`` is a :py:class:`bool` describing whether the representation contains complex conjugation (``True``) or not (``False``).

The following code creates the symmetries described above:

.. include :: ../../examples/TaAs2/symmetries.py
    :code: python
    :start-line: 8
    :end-line: 40


Expression basis
~~~~~~~~~~~~~~~~

Representation basis
~~~~~~~~~~~~~~~~~~~~



The :ref:`reference<reference>` gives you an overview of the available functions and classes.

Formalism
---------
Finally, let me give a more formal description of the problem. Let :math:`G` be the symmetry group that the Hamiltonian should respect, with a unitary representation :math:`D(g), g  \in G`. This imposes the symmetry constraint 

.. math ::
    \mathcal{H}(\mathbf{k}) = D(g) \mathcal{H}(g^{-1} \mathbf{k}) D(g^{-1}), \forall g \in G
    
on the :math:`\mathbf{k} \cdot \mathbf{p}` Hamiltonian. In the following, we will define a vector space containing :math:`\mathcal{H}`, and see how the symmetry constraint restricts the Hamiltonian to a certain subspace.

We want to consider only a certain form of :math:`\mathbf{k}` - dependence for the Hamiltonian, for example up to second order. So let :math:`V \subset \mathcal{F}(\mathbb{R}^3, \mathbb{R})` be the vector space which contains these functions of :math:`\mathbf{k}`. We require that :math:`V` is closed under

.. math ::
    \hat{F}_g: f \longmapsto \tilde{f}_g \\
    \tilde{f}(\mathbf{k}) = f(g^{-1}\mathbf{k}).

for all :math:`g \in G`. That is, 

.. math ::
    \forall g \in G, f \in V: \hat{F}_g(f) \in V.

:math:`\hat{F}_g` is a linear operator on :math:`V`. 

Let :math:`W` be the vector space of hermitian :math:`N \times N` matrices. 

.. math ::
    \hat{G}_g: ~&W & \longrightarrow W\\
    & A & \longmapsto D(g)A D(g^{-1})

is a linear operator on `W`. It is unitary under the Frobenius inner product.

Since the Hamiltonian is a hermitian matrix with :math:`\mathbf{k}`-dependence as given by :math:`V`, it follows that :math:`\mathcal{H} \in V\otimes W`. The symmetry constraints mean that 

.. math ::
    \forall g \in G:~ \left(\hat{F}_g \otimes \hat{G}_g\right)(\mathcal{H}) = \mathcal{H},

and thus

.. math ::
    \mathcal{H} \in \bigcap_{g\in G} \text{Eig}(\hat{F}_g \otimes \hat{G}_g , 1).

In conclusion, the problem of finding the general form of the Hamiltonian is equivalent to calculating this subspace.


.. toctree::
    :hidden:
    
    Usage and Formalism <self>
    example.rst
    reference.rst
