# Journal Club bot

Source code for the @JC_pathogenomic (https://twitter.com/JC_pathogenomic)

# Set up

**1. Clone the repo:**

    $ git clone https://github.com/happykhan/journalbot.git

**2. Install Python packages:**

    $ python -m pip install -r requirements.txt

**3. Set environment variables for Twitter authorisation (see https://apps.twitter.com/):**

Create a new twitter app, and generate a `consumer_key`, `consumer_secret`, `access_token`, `access_token_secret`. Export these as Environment Variables.
Also create an [NCBI/Pubmed Entrez email](https://www.ncbi.nlm.nih.gov/account/?back_url=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fgquery%2F).

Create a `.env` file in this directory, see `template.env`:

    consumer_key='xxx'
    consumer_secret='xxx'
    access_token='xxx'
    access_token_secret='xxx'
    entrez_email='user@email.com'

# Composing a search term

The bot will search for queries contained within the `config/seartchterms.txt` file. These should be composed in PubMed query syntax

- see https://www.ncbi.nlm.nih.gov/books/NBK3837/#_EntrezHelp_Entrez_Searching_Options_)

An example file would look like (use boolean operators and parentheses to evaluate more
complex searches):

     Smith J
     Smith J[author]
     Jones K
     E coli
     polymer AND synthesis
     nucleus NOT atomic
     human[organism] AND topoisomerase[protein name]
     Smith J OR Jones K
     Smith John M [FAU]

# Troubleshooting

Installation is easy with a Miniconda/Anaconda install:

    $ wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && bash Miniconda2-latest-Linux-x86_64.sh

Follow installation prompts. Allow appending to `$PATH` in `.bashrc`.

    $ which python # should yield: ~/anaconda2/bin/python or similar

should return a path including miniconda. If not, it's using the wrong python binary still, and you need to alter the order of paths in $PATH.
