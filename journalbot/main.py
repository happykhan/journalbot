#! /usr/bin/env python
"""
Simple Journal club twitter bot

It picks a set of user defined authors and then proceeds to post twitter
updates.

### CHANGE LOG ###
2014-12-23 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Initial build
2014-12-28 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Added basic functionality. Scirpt loads papers from text files,
      updates to twitter, fetches paper metadata.
2016-11-08 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Patched libraries, timeout handling
    * Replaced citation counters with IF factoring
2017-07-20 Joe Healey <jrj.healey@gmail.com>
    * Update biopython requirement to 1.70
    * Use of credentials file
    * Bug fixes for pubmed search
    * Remove specificity for authorship to allow
      full PubMed query syntax including keywords
2018-04-09 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Modified paper selection system. Papers now
      must be within the 30 days. Ranked by date
      & IF
    * Cleans up random html tags in titles.
    * Post with author as twitter handle. Changes to
      authorlist file!
2022-11-14 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Major refactor
"""
import datetime
import csv
import tweepy
import os
import logging
import time
import re
from Bio import Entrez
from journalbot.util import generate_message, statushash, calculateScore, cleanhtml, load_credentials
import pytz

def get_next_paper(paperlist, journal_file, search_file, CUTOFF=60):
    Entrez.email = os.environ.get("entrez_email", "nabil@happykhan.com")
    # Retrieve all papers for each author from file
    # Use entrez pubmed search
    logging.info("Fetching next paper to post")
    # Retrieve file saved last update time if available
    search_terms = {}
    total_string = ""
    # Handle google sheet as well.
    term_list = []
    journal = {}
    ignore_terms = [] 
    for term in open('config/ignore_terms.txt').readlines():
        ignore_terms.append(term.strip())
    for term in csv.DictReader(open(search_file), delimiter=";"):
        term_list.append(term["search_term"].strip())        
    for journal_row in csv.DictReader(open(journal_file), delimiter="\t"):
        journal[journal_row['name']] = journal_row['count']
    for term in term_list:
        total_string += f" OR ({term})"
        position = round(len(total_string) / 2500)
        if search_terms.get(position):
            search_terms[position] += f" OR ({term})"
        else:
            search_terms[position] = f"({term})"
    entries = []
    for search_string in search_terms.values():
        try:
            logging.info("Fetching papers relating to " + search_string)
            handle = Entrez.esearch(db="pubmed", term=search_string, reldate=CUTOFF, retmax=4000)
            record = Entrez.read(handle)
            handle.close()
            time.sleep(3)
            if len(record["IdList"]) > 0:
                entries += record["IdList"]
        except Exception:
            logging.error("Failed to fetch " + search_string + ". Skipping.")
    chunks = [entries[x : x + 1000] for x in range(0, len(entries), 1000)]
    for entrylist in chunks:
        result = []
        try:
            handle = Entrez.esummary(db="pubmed", id=",".join(entrylist))
            result = Entrez.read(handle)
            handle.close()
            for r in result:
                try:
                    pub_count = journal.get(r["FullJournalName"].lower(), 0)
                    ignore = False
                    for x in ignore_terms:
                        if re.search(x, r["Title"]):
                            ignore = True

                    if not ignore:
                        newpaper = dict(
                            authors=r["AuthorList"],
                            pmid=r["ArticleIds"]["pubmed"][0]
                            .encode("ascii", "ignore")
                            .decode("utf-8"),
                            full_title=cleanhtml(r["Title"].encode("ascii", "ignore")),
                            journal=r["FullJournalName"],
                            pub_count=pub_count,
                            date=r["EPubDate"],
                        )
                        twit_length = 0
                        newpaper["tweet_title"] = newpaper["full_title"][
                            0 : (280 - (twit_length + 20))
                        ]
                        newpaper["score"] = calculateScore(newpaper, CUTOFF)
                        newpaper["status_hash"] = statushash(newpaper["tweet_title"])
                        # Is the paper already there?
                        if (
                            newpaper["score"] > 0
                            and newpaper.get("pmid")
                            and not newpaper["status_hash"]
                            in (item["status_hash"] for item in paperlist)
                        ):
                            paperlist.append(newpaper)
                except:
                    logging.error(
                        "Error in loading %s" % r["Title"].encode("ascii", "ignore")
                    )
        except:
            logging.error("Error in loading %s" % entrylist[0])
    paperlist.sort(key=lambda paper: paper["score"], reverse=True)
    paperlist_path = os.path.join(os.path.dirname(journal_file), 'paperlist.tsv')
    out_file = csv.DictWriter(
        open(paperlist_path, "w", newline=""),
        fieldnames=['pmid', 'full_title', 'journal', 'pub_count', 'date', 'tweet_title', 'status_hash', 'posted', 'score'],
        delimiter="\t",
    )    
    out_file.writeheader()
    for paper in paperlist: 
        new_row = paper.copy()
        new_row.pop('authors', None)
        out_file.writerow(new_row)
    logging.info("Written papelist file to %s" %paperlist_path)
    for next_paper in paperlist:
        if next_paper["score"] > 0:
            logging.info(generate_message(next_paper))
    return paperlist[0]


def find_journal(search_file, journal_file):
    Entrez.email = os.environ.get("entrez_email", "nabil@happykhan.com")
    journals = {}
    search_terms = {}
    total_string = ""
    term_list = []
    for term in csv.DictReader(open(search_file), delimiter=";"):
        term_list.append(term["search_term"].strip())
    for term in term_list:
        total_string += f" OR ({term})"
        position = round(len(total_string) / 2500)
        if search_terms.get(position):
            search_terms[position] += f" OR ({term})"
        else:
            search_terms[position] = f"({term})"
    entries = []
    for search_string in search_terms.values():
        try:
            logging.info("Fetching papers relating to " + search_string)
            handle = Entrez.esearch(db="pubmed", term=search_string, retmax=4000)
            record = Entrez.read(handle)
            handle.close()
            time.sleep(3)
            if len(record["IdList"]) > 0:
                entries += record["IdList"]
        except Exception:
            logging.error(f"Failed to fetch {search_string}. Skipping.")
    chunks = [entries[x : x + 1000] for x in range(0, len(entries), 1000)]
    for entrylist in chunks:
        result = []
        try:
            handle = Entrez.esummary(db="pubmed", id=",".join(entrylist))
            result = Entrez.read(handle)
            handle.close()
            for r in result:
                journal_name = r["FullJournalName"].lower()
                if journals.get(journal_name):
                    journals[journal_name]["count"] += 1
                else:
                    journals[journal_name] = dict(count=1, name=journal_name)
        except:
            logging.error("Error in loading %s" % entrylist[0])
    out_file = csv.DictWriter(
        open(journal_file, "w", newline=""),
        fieldnames=["name", "count"],
        delimiter="\t",
    )
    out_file.writeheader()
    out_file.writerows(journals.values())
    logging.info("Written journal file to %s" %journal_file)
    return journals


def journal_thread(search_file, journal_file):
    logging.info("Starting journal thread")
    while True:
        logging.info("Running journal update")
        find_journal(search_file, journal_file)
        time.sleep(86400)  # Wait one day


def post_thread(search_file, journal_file, tinterval, date_cutoff):
    logging.info("Running posting thread")
    (
        consumer_key,
        consumer_secret,
        access_token,
        access_token_secret,
    ) = load_credentials()
    # OAuth process, using the keys and tokens
    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
    auth.set_access_token(access_token, access_token_secret)

    # Creation of the actual interface, using authentication
    api = tweepy.API(auth)
    lastTweetTime = datetime.datetime.utcfromtimestamp(0)
    for status in tweepy.Cursor(api.user_timeline).items(1):
        lastTweetTime = status.created_at

    # Retrieve existing statuses
    while True:
        logging.info("Waiting. Last Tweet at " + str(lastTweetTime))
        if (
            pytz.utc.localize(datetime.datetime.now())
            > (lastTweetTime + datetime.timedelta(minutes=tinterval))
            and os.path.exists(journal_file)
        ):
            logging.info("Fetching previous tweets")
            page_list = []
            for page in tweepy.Cursor(api.user_timeline, count=200).pages(10):
                page_list.append(page)
            # Retrieve previous tweets
            title_match = re.compile("(NEW: |OF NOTE: )?(.+) http.+")
            paperlist = []
            for page in page_list:
                for status in page:
                    # Tweepy seems to truncate tweets by default.
                    if status.truncated:
                        clean_status = (
                            api.get_status(status.id, tweet_mode="extended")
                            ._json["full_text"]
                            .encode("ascii", "ignore")
                        )
                    else:
                        clean_status = status.text.encode("ascii", "ignore")
                    clean_status = clean_status.decode("utf-8")

                    if title_match.match(clean_status):
                        tweet_title = title_match.match(clean_status).group(2)
                        status_hash = statushash(tweet_title)
                        if not status_hash in (
                            item["status_hash"] for item in paperlist
                        ):
                            paperlist.append(
                                dict(
                                    tweet_title=tweet_title,
                                    score=0,
                                    posted=True,
                                    status_hash=status_hash,
                                )
                            )
            logging.info("Updated previous tweets")
            logging.info("Last Tweet at " + str(status.created_at))
            next_paper = get_next_paper(
                paperlist,
                journal_file,
                search_file,
                date_cutoff,
            )
            if next_paper["score"] > 0:
                try:
                    message = generate_message(next_paper)
                    # api.update_status(message)
                    # lastTweetTime = datetime.datetime.now()
                except:
                    logging.error("Error posting status: %s" % message)
            else:
                # All possible posts have score 0, Nothing new to post!
                lastTweetTime = datetime.datetime.now()
        time.sleep(600)
