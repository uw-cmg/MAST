###################################
MAST development workflow
###################################

Last modified TTM 2016-02-11

***********************************
I am a developer
***********************************

=======================
Create a new branch
=======================

Go to the repository at `<http://git@github.com/uw-cmg/MAST>`_

Create a new branch from dev.

.. _create-test-env:

=========================================================   
Create new anaconda test environment with dependencies
=========================================================
In the instructions below, the test environment is called ``mast_myname_test``.::

    conda create -n mast_myname_test numpy scipy matplotlib pandas nose sphinx
    source activate mast_myname_test
    
    (mast_myname_test)[user@cluster ~]$ which conda
    ~/anaconda/envs/mast_myname_test/bin/conda
    (mast_myname_test)[user@cluster ~]$ which python
    ~/anaconda/envs/mast_myname_test/bin/python
    (mast_myname_test)[user@cluster ~]$ which pip
    ~/anaconda/envs/mast_myname_test/bin/pip

Install MAST dependencies by pip-installing MAST::
    
    (mast_myname_test)[user@cluster ~]$ pip install --cert=~/anaconda/ssl/cacert.pem MAST

The cert tag may be important for finding newer versions etc.

If unsure, use pip -v -v -v (v for verbose)

You may also have to manually install ase and add it to PYTHONPATH

Now, you would actually want to develop on a git clone version of MAST rather than on the package-installed MAST.

So, unlink the package-installed MAST::

    cd ~/anaconda/envs/mast_myname_test/lib/python2.7/site-packages/
    mv MAST MAST-do-not-find

Clone a fresh copy and link to your own branch of MAST::

    cd
    git clone ssh://git@github.com/uw-cmg/MAST --branch myname mast_myname

The instructions above would put your MAST repository clone, with the branch name ``myname``, into the folder ``//home/<user>/mast_myname``

Add your repository to the activation script::

    cd ~/anaconda/envs/mast_myname_test/bin
    vi activate

...using a line like this::

    export PYTHONPATH=$HOME/mast_myname:$PYTHONPATH

Make sure the repository is getting your branch version::

    source deactivate
    source activate mast_myname_test
    python
    >>> import MAST
    >>> help(MAST)

Look for ``FILE`` to be ``//home/<user>/mast_myname/MAST/__init__.py``

Run unit tests to check the installation. Occasionally nosetests will use the wrong python, so if that happens, instead of using the command ``nosetests``, use::

    cd //home/<user>/mast_myname/MAST/test
    python -m nose

or::

    python -m nose -v --nocapture

========================
Code on your branch
========================

Code on your branch in ``//home/<user>/mast_myname``
Create your own unit tests in ``//home/<user>/mast_myname/MAST/test``

Commit and push changes to your branch.::

    git commit -a --message="commit message"
    git push ssh://git@github.com/uw-cmg/MAST myname:myname

Run unit tests by::

    source ~/anaconda/bin/activate mast_myname_test
    cd ~/mast_myname/MAST/test
    python -m nose

(Modify ``test/workflow_test/WORKFLOW_CONFIG`` and ``test/workflow_test/examples`` in order to run full workflow tests.)

You can also run unit tests singly, running only those tests directly affected by your changes. See nosetests documentation at python.org.

After development, run all unit tests.

Investigate failed tests. Repair code and/or tests as necessary.

When all tests have passed, commit and push to own branch.

.. _reconcile-with-dev:

===================================================   
Reconcile your branch with the development branch
===================================================

Pull from dev to your local clone::

    git pull ssh://git@github.com/uw-cmg/MAST dev

Address any merge errors.

Run unit tests.

Investigate failed tests. Repair code and/or tests as necessary, consulting with other developers as necessary.

Commit and push to your branch on github.

When all tests have passed, commit and push to dev.::

    git push ssh://git@github.com/uw-cmg/MAST myname:dev

.. _release-manager:

*********************************
I am the release manager
*********************************

===============================
Create a dev-test environment
===============================

Repeat the steps above under :ref:`create-test-env` with a new environment and the dev branch.

For the following instructions, assume your new environment is ``mast_dev_test``, and that you cloned ``--branch dev`` into the folder ``mast_dev``.

If any unit tests fail at this stage, assign any problems to developers.
When developers have finished fixing tests, pull from dev again.

Each developer starts over from the :ref:`reconcile-with-dev` step.

=========================================
Check documentation and upload to pypi
=========================================

Change the version number in ``//home/<user>/mast_dev/MAST/_version.py``

Edit the sphinx-build path in ``//home/<user>/mast_dev/MAST/doc/Makefile`` to correspond to ``//home/<user>/anaconda/envs/mast_dev_test/bin/sphinx-build``.

Build the documentation::

    source ~/anaconda/bin/activate mast_dev_test
    (mast_dev_test)[user@cluster ]$ cd //home/<user>/mast_dev/MAST/doc
    (mast_dev_test)[user@cluster ]$ make html

Investigate any bad builds and edit the rst files in ``/MAST/doc/source`` as necessary.

Copy the MAST/doc/build/html folder to your desktop.

Preview the html.

Modify ``MAST/doc/source`` pages as appropriate.

Rebuild and preview until the documentation is complete.

Zip the entire contents of html folder into a zip file by selecting all the files and folders in the html folder and compressing them (this has to be done from the desktop for some reason, or else it will not be readable by pythonhosted.org)

Log in to the pypi page at https://pypi.python.org/pypi/MAST
At the bottom of the page, upload the zip into the Upload documentation section
The upload may spin for a while, but it should work.

============================================================
Package to the pypi test repository and test the package
============================================================

Create a ``~/.pypirc`` file that looks something like::

    [distutils]
    index-servers=
        pypi
        test

    [pypitest]
    repository = https://testpypi.python.org/pypi
    username = <username>
    password =  <pwd>

    [pypi]
    repository = https://pypi.python.org/pypi
    username = <username>
    password =  <pwd>

Go to ``//home/<user>/mast_dev``
Change ``MAST/_version.py`` to a use a dummy number (e.g. if the current version is 1.2.1, use 1.2.0.1)
Upload to the test pypi::

    cd //home/<user>/mast_dev
    python setup.py register -r pypitest sdist upload -r pypitest

If you are in a conda environment, deactivate it::
    
    source deactivate
    
Make a new conda environment::

    conda create -n mast_testpypi_test numpy scipy matplotlib pandas nose sphinx
    source activate mast_testpypi_test

**If there is** a line adding a specific version of MAST in your ~/anaconda/envs/<new environment>/activate script, remove it, LOG OUT, and source activate your new environment. (Otherwise you will not install MAST from the pypi test repository.)

To install MAST from the test server::

    pip install --cert ~/anaconda/ssl/cacert.pem pymatgen
    pip install --cert ~/anaconda/ssl/cacert.pem custodian
    pip install -i https://testpypi.python.org/pypi --cert ~/anaconda/ssl/cacert.pem MAST

(It is necessary to install some dependencies separately if they are not on the test server.)

Do not unlink MAST as you did for the mast_dev_test environment. This time, you want the package-installed MAST version.

Run unit tests.

Eventually this should work::
    
    (mast_testpypi_test)[user@cluster ~]$ cd ~/anaconda/envs/mast_testpypi_test/lib/python2.7/site-packages/MAST/test
    (mast_testpypi_test)[user@cluster ~]$ python -m nose

But in practice not all files get copied over, so use the test files from dev (which should be up to date, since you created the testpypi package from dev) but use this packaged version's python and MAST.::

    cd //home/<user>/mast_dev/mast_20160211_dev/MAST/test
    python -m nose

If any unit tests fail, go back to the developers.

Each developer starts over from the :ref:`reconcile-with-dev` step.

When all unit tests have passed, the release manager starts over from :ref:`release-manager`.

If all tests pass, package to real repository.

========================================
Package to the real pypi repository
========================================
Change ``~/mast_dev/MAST/_version.py`` to the correct version.

Package to the real pypi repository::

    cd //home/<user>/mast_dev
    python setup.py register -r pypi sdist upload -r pypi

Make a new environment as above and source it::

    source deactivate
    conda create -n mast_pypi_test numpy scipy matplotlib pandas nose sphinx
    source activate mast_pypi_test
    pip install --cert ~/anaconda/ssl/cacert.pem MAST

Run unit tests out of ``~/mast_dev``.

Go back to developers if any unit tests fail. Nothing should fail from here, but it might.

If all tests pass, then only archiving, cleanup, and notifications are left.

============
Archiving
============

Sign in to github.com

* Draft a new release
* Use the same version number as in the pypi real repository as a tag (conventionally, v#.#.#)
* Use MAST-v#.#.# as the title
* Write a small description
* Publish release
* Save the .tar.gz file for this release

Sign in to zenodo.org

* Create a new Software upload with the .tar.gz file

=====================
Cleanup
=====================

* Remove all created environments using ``conda env remove -n <environment name>``
* Remove all cloned MAST directories (except any that you use for personal coding)
* Remove any unneeded github branches that may have been created during the testing phase (e.g. if two developers create a mutual branch in order to test some code conflict)

=================
Notification
=================
Send out an email to the MAST development team and any core users. Done!

