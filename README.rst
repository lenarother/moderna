================================================
ModeRNA - A program for comparative RNA modeling
================================================

Copyright 2009 by Magdalena Rother, Kristian Rother, Tomasz Puton and Janusz M. Bujnicki

Version 1.7.1

Homepage:
   iimcb.genesilico.pl/moderna

Technical Support:
   rother.magdalena@gmail.com



Installation Instructions
=========================

0. Quick guide
--------------
 
   python setup.py install

1. Requirements
---------------

ModeRNA runs on any modern Windows or Linux PC. There are two versions available:

- Source code - full functionality, but requires some libraries.
- Windows Binary - easy to install, but limited functionality.


2. Installing the Windows Binary
-------------------------------- 
All you need to do is:

- Download the latest binary [Moderna Version 1.5.0 (.exe)].

- Unzip the archive.

- Open a command line window and switch to the unzipped directory. Type 

  python moderna.py -h 

  to see available options.


3. Installation on Linux
------------------------

To install ModeRNA on Linux, you need to:

- Download the source distribution [Moderna Version 1.5.2 (source)].

- Unzip the archive. A catalog with the main program moderna.py is created.

- Make sure Python 2.5 or a higher version is installed. (on Ubuntu Linux, use sudo apt-get install python.

- Make sure Numpy is installed (sudo apt-get install python-numpy).

- Make sure BioPython is installed (sudo apt-get install python-biopython).
   (tested with BioPython 1.51-1.53)

- Make sure the moderna/moderna.py file is executable:
  chmod a+x moderna/moderna.py

- Add the path to the moderna directory to your PYTHONPATH variable,e.g. 

- run:
  python setup.py install

- After this, you can from any location write:

python
>>> from moderna.moderna import *


4. Installing ModeRNA on Windows
--------------------------------

To install ModeRNA on Windows, you need to:

- Download the source distribution [Moderna Version 1.5.2 (source)].

- Unzip the archive. A catalog with the main program moderna.py is created.

- Make sure that these libraries are installed:
    (all available from iimcb.genesilico.pl/moderna)
    Python 2.5 or a higher
    Numpy
    BioPython (tested with BioPython 1.51-1.53)

- Run from the terminal:
  C:/Python25/python.exe setup.py install

- After this, you can write in the Python shell:

>>> from moderna.moderna import *


Legal Disclaimer
----------------

ModeRNA is released under the GPL license, a copy of which is included in 
the distribution (See LICENSE_GPL.TXT for details). For the files in the 
PDB/ directory, the Biopython License applies as well. 
See PDB/LICENSE_BIOPYTHON.TXT for details).

This software is provided "as-is". There are no expressed or implied 
warranties of any kind, including, but not limited to, the warranties of 
merchantability and fitness for a given application. In no event shall 
the authors be liable for any direct, indirect, incidental, special, 
exemplary or consequential damages (including, but not limited to, loss 
of use, data or profits, or business interruption) however caused and on 
any theory of liability, whether in contract, strict liability or tort 
(including negligence or otherwise) arising in any way out of the use 
of this software, even if advised of the possibility of such damage.

The authors take no responsibility for damage caused by this program 
or its components. 


Contributors
------------

- Magdalena Rother   - implementation
- Pawel Piatkowski   - implementation
- Kristian Rother    - architecture and unit tests
- Tomasz Puton       - model validation and testing
- Janusz Bujnicki    - concept and supervision


Acknowledgements
----------------

Credit goes to our lab colleagues Pawel Skiba, Piotr Byzia, Irina Tuszynska, 
Joanna Kasprzak, Jurek Orlowski, Pawel Lukasz, Tomasz Osinski, Marcin 
Domagalski, Anna Czerwoniec, Stanislaw Dunin-Horkavic, Marcin Skorupski, 
and Marcin Feder for their comments and constructive criticism during 
development. 

The PDB parser ued by Moderna uses BioPython with kind support by 
Thomas Hamelryck. The unit test framework was brought near to us by 
Sandra Smit, Rob Knight, and Gavin Huttley. We also would like to thank 
Neocles Leontis, Fabrice Jossinet, Francois Major, and Eric Westhof who 
provided helpful advice on various occasions.
Special thanks go to the group of Russ Altman, who provided us with 
their modeling example to test ModeRNA.


References
----------

Components of ModeRNA are based upon the following pieces of scientific literature:

- [1] Czerwoniec A, Dunin-Horkawicz S, Purta E, Kaminska KH, Kasprzak JM, Bujnicki JM, Grosjean H, Rother K. MODOMICS: a database of RNA modification pathways. 2008 update. Nucleic Acids Res. 2008 Oct 14.
- [2] Yang H, Jossinet F, Leontis N, Chen L, Westbrook J, Berman H, Westhof E. Tools for the automatic identification and classification of RNA base pairs. Nucleic Acids Res. 2003 Jul 1;31(13):3450-60.
- [3] Gendron P, Lemieux S, Major F. Quantitative analysis of nucleic acid three-dimensional structures. J Mol Biol. 2001 May 18;308(5):919-36.
- [4] Michalsky E, Goede A, Preissner R. Loops In Proteins (LIP)–a comprehensive loop database for homology modelling. Protein Eng. 2003 Dec;16(12):979-85.
- [5] Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422-3. Epub 2009 Mar 20.

-------------------------------------------------------------------------


============================
Notes from pyscaffold README
============================

This project was set up with PyScaffold 0.9.
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


