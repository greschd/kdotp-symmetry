.. _home:

kdotp-symmetry
==============

This is a tool to calculate the general form of a :math:`\mathbf{k}\cdot\mathbf{p}` Hamiltonian under a given symmetry constraint.

Usage
-----

You can install this tool with with pip:

.. code ::

    pip install kdotp-symmetry

Its usage is best explained with an example. 


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
    
    Usage and Formalism <self>
    reference.rst
