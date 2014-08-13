=======
moderna
=======

Your project was successfully set up with PyScaffold 0.9.
Following features are available:

Packaging
=========

Run ``python setup.py sdist``, ``python setup.py bdist`` or
``python setup.py bdist_wheel`` to build a source, binary or wheel
distribution.


Complete Git Integration
========================

Your project is already an initialised Git repository and ``setup.py`` uses
the information of tags to infer the version of your project with the help of
`versioneer <https://github.com/warner/python-versioneer>`_.
To use this feature you need to tag with the format ``vMAJOR.MINOR[.REVISION]``
, e.g. ``v0.0.1`` or ``v0.1``. The prefix ``v`` is needed!
Run ``python setup.py version`` to retrieve the current `PEP440
<http://www.python.org/dev/peps/pep-0440/>`_-compliant version. This version
will be used when building a package and is also accessible through
``my_project.__version__``.
The version will be ``unknown`` until you have added a first tag.


Sphinx Documentation
====================

Build the documentation with ``python setup.py docs`` and run doctests with
``python setup.py doctest``. Start editing the file ``docs/index.rst`` to
extend the documentation.


Unittest & Coverage
===================

Run ``python setup.py test`` to run all unittests defined in the subfolder
``tests`` with the help of `py.test <http://pytest.org/>`_. The py.test plugin
`pytest-cov <https://github.com/schlamar/pytest-cov>`_ is used to automatically
generate a coverage report. For usage with a continuous integration software
JUnit and Coverage XML output can be activated. Checkout ``putup -h`` for
details.

Requirements Management
=======================

Add the requirements of your project to the ``requirements.txt`` file which
will be automatically used by ``setup.py``.


Easy Updating
=============

Keep your project's scaffold up-to-date by applying
``putput --update my_project`` when a new version of PyScaffold was released.
It may also be used to change the url, license and description setting.


.. note::

    Replace the content of this file with a description of your project.

