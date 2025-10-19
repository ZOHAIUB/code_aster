.. _devguide-naming-conventions:

*******************************
Naming conventions / Python API
*******************************

.. note::
    This section is still a *work in progress* document.


=============
General rules
=============

A general rule is to follow the `CamelCase <https://en.wikipedia.org/wiki/Camel_case>`_
convention in C++ code *and* its Python Bindings.
By extension, the *high level* user API (Python) also follows *CamelCase* naming.
The rest of the Python code, *low level* functions, should follow the
`PEP8 Style Guide for Python code <https://www.python.org/dev/peps/pep-0008/>`_
that does not recommend *CamelCase* naming.

- The names of **classes** start with an upper case (example: *MeshCoordinatesField*). The shared_ptr on this class has the same name with Ptr at the end (example: *MeshCoordinatesFieldPtr*)

- The names of **methods**/**members** start with a lower case and
  does not contain the object name (example: *getValues* and not *getFieldValues*).

- The names of specialized **methods**/**classes** depending on the type have to add this type always at the end. (exemple: *FieldOnNodesReal* and *FieldOnNodesComplex* for classes, and *getValueReal* and *getValueVectorOfReal* for methods)

- Reuse existing names in newly created objects. Have a look to the :ref:`genindex` page
  should give good ideas for new methods.
  One can also search for names in the :ref:`devguide-objects_Datastructure` page.

- Plural forms should be used when it is relevant.

- Use same method name to pass a single value or a vector.
  Example: Do not define :py:obj:`addMaterialOnMesh(mater)` and
  :py:obj:`addMaterialsOnMesh(vect_mater)`, but only
  :py:meth:`~code_aster.Objects.MaterialField.addMaterialsOnMesh` with the both
  interfaces.

- Use same method name to return a single value or a vector if the argumets are differents.
  Prefer singular version.
  Example: Do not define :py:obj:`getNodesFromDOFs()` and
  :py:obj:`getNodeFromDOF(DOFId)`, but only
  :py:meth:`~code_aster.Objects.EquationNumbering.getNodeFromDOF` with the both
  interfaces.

- Use same method name to pass a differents values.
  Example: Do not define :py:meth:`addRealValue( realVal )` and
  :py:meth:`addComplexValue(cmplxVal)`, but only
  :py:meth:`addValue` with the both
  interfaces.

- At least at the beginning, only pure Python objects are returned (example: *list* or
  *list[list]* and not *numpy* arrays).

- *std::vector* and *JeveuxVector* objects are automatically converted to Python *list*.
  *JeveuxCollection* objects are also converted as *list[list]*.
  See :file:`PythonBindings/ConvertersInterface` for supported converters.

- Strings values are returned without trailing spaces.

- build : method to construct or modify internal objects of the class. This method return *void* or *bool*

- updateValuePointers : method to update all JeveuxVector of a class

- computeSomething : a method begins with compute if it return an object which is explicitely computed and not an internal object (exemple: computeDirichletBC, computeLoads)

- prefer FromSomethong to AssociatedToSomething: this is shorter and explicit.

Another rule is to define elementary methods not to create objects with a huge size.
For example, the groups names can be extracted from a mesh and cells indexes can be
extracted for a named group. But no method returns a *dict* of cells indexes for all
groups for example.
Another example: do not create shortcuts as *mesh.getCoordinatesValues()*
but *mesh.getCoordinates().getValues()*.

A senior developer must validate any changes to the user API (C++/Python bindings and
pure Python API). Do not hesitate to ask a precion to a senior developper for C++/Python and complete this documentation


========
Glossary
========

Mesh object
-----------

Terms for the :py:class:`~code_aster.Objects.Mesh` object:

- *mesh* is an object composed of *nodes* and *cells*.
- *node* for a node.
- *cell* for an element of the mesh (*element* term is used for a *finite element*). It was "Maille" in old french code_aster terminology
- *virtual nodes* for a node not originally present in the mesh (it is used to create loads). It was "Noeud tardif" in old french code_aster terminology
- *virtual cell* for a cell not originally present in the mesh (it is used to connect nodes or virtual nodes). It was "Elément tardif" in old french code_aster terminology
- *GroupOfNodes* for a group of nodes, named with at most 24 chars.  It was "GROUP_NO" in old french code_aster terminology
- *GroupOfCells* for a group of cells, named with at most 24 chars.  It was "GROUP_MA" in old french code_aster terminology
- *Connectivity* for the mesh connectivity.
- For a *ParallelMesh*, an additional boolean argument named *local* allows to work
  on the local part (that belongs to each MPI process, *local=True*) or on the
  global mesh (*local=False*).

Methods are applied on all the mesh: *OnMesh*, on a group of cells: *OnGroupOfCells*
or on a group of nodes *OnGroupOfNodes*.

.. todo::
    Add *same* methods to *ParallelMesh* with a *local* argument.


Model object
------------

- *element* for a finite element (not a *cell*).


Result objects
--------------

- *result* is an object that contains several fields and eventually some other properties.


Numbering objects
-----------------

- *DOF* is a degree of freedom (*DOFs* for the plural form").
- *LagrangeDOF* is a Lagrange DOF on a virtual node.