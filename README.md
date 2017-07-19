GibsonGroupBot
==============

Source code for the GibsonGroupBot @ https://twitter.com/GibsonGroupBot

Forked from @happykhan/journalbot

# Set up:

Breaks under python3 so ensure it's run with Python2.7.

1. Clone the repo:
------------------

    $ git clone https://github.com/jrjhealey/journalbot.git


2. Make sure all dependencies are satisfied:
-------------------------------------------

    $ python -m pip install -r /path/to/journalbot/requirements.txt


3. Set environment variables for Twitter authorisation (see https://apps.twitter.com/):
---------------------------------------------------------------------------------------
Create a new twiiter app, and generate a `consumer_key`, `consumer_secret`, `access_token`, `access_token_secret`. Export these as Environment Variables.
Also create an NCBI/Pubmed Entrez email (https://www.ncbi.nlm.nih.gov/account/?back_url=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fgquery%2F).

E.g. in a script:

    export consumer_key='xxxxxxxxxxxxxxxxxxxx'
    export consumer_secret='xxxxxxxxxxxxxxxxxxxx'
    export access_token='xxxxxxxxxxxxxxxxxxxxx'
    export access_token_secret='xxxxxxxxxxxxxxxxxxx'
    export entrez_email='user@email.com'

Then source this file to create the variables:

    $ source /path/to/credentials_file

# Troubleshooting

Installation is easy with a Miniconda/Anaconda install:

    $ wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && bash Miniconda2-latest-Linux-x86_64.sh

Follow installation prompts. Allow appending to `$PATH` in `.bashrc`.

    $ which python # should yeild: ~/anaconda2/bin/python or similar

should return a path including miniconda. If not, it's using the wrong python binary still, and you need to alter the order of paths in $PATH.
