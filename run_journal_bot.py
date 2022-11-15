#! /usr/bin/env python
"""
Simple Journal club twitter bot

It picks a set of user defined authors and then proceeds to post twitter
updates.

CLI starter script for journalbot
"""

import sys
import os
import traceback
import argparse
import time
import threading
import logging
from dotenv import load_dotenv
import journalbot.__init__ as meta
from journalbot.main import post_thread, journal_thread
from journalbot.util import load_credentials

epi = (
    "Licence: "
    + meta.__licence__
    + " by "
    + meta.__author__
    + " <"
    + meta.__author_email__
    + ">"
)


def main():

    global args
    # Init paper db update thread
    log_level = logging.INFO
    if args.verbose:
        log_level = logging.DEBUG
    logging.basicConfig(
        level=log_level, format="[%(levelname)s] (%(threadName)-10s) %(message)s"
    )
    threads = []
    logging.info("Starting JournalBot")
    search_file = os.path.join("config", "searchterms.txt")
    journal_file = os.path.join("config", "journal.tsv")
    if not os.path.exists(search_file):
        logging.error(
            "Search terms file not found at %s" %search_file
        )
        sys.exit(-1)

    # Check credentials
    load_dotenv()
    (
        consumer_key,
        consumer_secret,
        access_token,
        access_token_secret,
    ) = load_credentials()
    if not (consumer_secret and consumer_key and access_token and access_token_secret):
        logging.error(
            "TOKENS NOT FOUND - Populate the .env file."
        )
        sys.exit(-1)

    # Init Posting thread
    post = threading.Thread(
        target=post_thread,
        args=(
            search_file,
            journal_file,
            args.tinterval,
            args.date_cutoff,
        ),
    )
    threads.append(post)

    # Init Journal thread
    post = threading.Thread(
        target=journal_thread,
        args=(
            search_file,
            journal_file,
        ),
    )
    threads.append(post)
    for t in threads:
        t.start()
    for t in threads:
        t.join()


if __name__ == "__main__":
    try:
        start_time = time.time()
        desc = __doc__.split("\n\n")[1].strip()
        parser = argparse.ArgumentParser(description=desc, epilog=epi)
        parser.add_argument(
            "-v", "--verbose", action="store_true", default=False, help="verbose output"
        )
        parser.add_argument(
            "--version", action="version", version="%(prog)s " + meta.__version__
        )
        parser.add_argument(
            "-t",
            "--tinterval",
            action="store",
            help="Time interval between tweets (minutes)",
            type=int,
            default=300,
        )
        parser.add_argument(
            "--workdir", action="store", help="Working directory", default="config"
        )
        parser.add_argument(
            "-d",
            "--date_cutoff",
            action="store",
            help="Post papers in the last X days [def: 60]",
            type=int,
            default=60,
        )
        args = parser.parse_args()
        if args.verbose:
            print("Executing @ " + time.asctime())
        main()
        if args.verbose:
            print("Ended @ " + time.asctime())
        if args.verbose:
            print("total time in minutes: " + time.time() - start_time / 60.0)
        sys.exit(0)
    except Exception:
        print("ERROR, UNEXPECTED EXCEPTION")
        traceback.print_exc()
        os._exit(1)
