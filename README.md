optmage
=======

Scripts for designing oligos and primers for MAGE (Multiplex Automated Genome Engineering)

## Requirements

The following are third-party software packages that must be installed on the system. Python requirements are described in the Installation section below.

### Unafold (required for secondary structure calculation).

1. Download and install from here
   [http://dinamelt.rit.albany.edu/download.php](http://dinamelt.rit.albany.edu/download.php)

* If using linux, you probably want to download the RPM and use these instructions:
<https://overlappingminds.com/sh/thoughts/0aa6d79e-fb8b-4f84-b287-f8e4494eac49>


## Installation

The setup.py config isn't ready yet, so please use the following installation procedure:

1. Install [virtualenv](http://www.virtualenv.org/en/latest/index.html) if you don't have it yet. (You may want to install [pip](http://pypi.python.org/pypi/pip/) first.)

2. Create a new virtual environment for this project. This virtual environment isn't part of the project so just put it somewhere on your machine.

        $ virtualenv /path/to/venv/project-venv

    If you want to use a version of python different from the OS default you can specify the python binary with the '-p' option:

        $ virtualenv -p /usr/local/bin/python2.7 /path/to/venv/project-venv

3. Activate the environment in the shell. This will use `python` and other binaries like `pip` that are located your pyenv. You should do this whenever running any python/django scripts.

        $ source /path/to/venv/project-venv/bin/activate .

4. Install the dependencies in your virtual environment. We've exported the requirements in the requirements.txt file. In theory, these should all be installable with the single command:

        (venv)$ pip install -r requirements.txt

However, in reality, this doesn't seem to work perfectly. Specific issues:

* Manually install numpy using pip first.

        (venv)$ pip install numpy

## Running Tests

With the correct python environment, run:

    $ python setup.py test

If the above isn't working, you can use `nose` directly for testing from the top-level directory (one directory **above** the `tests/` directory), by running (with the virtualenv enabled on the shell):

    (venv)$ nosetests


## Usage

There is currently a single script that does all the work: `src/optmage/oligo_designer.py`.

The user can list all the options available by running (remember to use the proper virtualenv):

    (venv)$ python oligo_designer.py -h

By default the script creates oligos based on the configuration in `src/optmage/data/INPUT_targets.txt`, relative to the reference genome `src/optmage/data/mg1655.fasta`. Users should specify different inpus using the above flags.

NOTE: The script currently has defaults tuned to the canonical E. coli MG1655 strain so using a different strain may require some edits.
