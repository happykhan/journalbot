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
"""

import sys, os, traceback, argparse, warnings
import time, re
import __init__ as meta 
import threading
import logging
from time import gmtime
import datetime
import csv
import requests
import json
import tweepy
from Bio import Entrez
import gspread

epi = "Licence: "+meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"

ifactor = {} 
blacklist = {} 

consumer_key = None
consumer_secret = None
access_token =  None
access_token_secret = None

def main ():

    global args
    # Init paper db update thread 
    log_level = logging.INFO
    if args.verbose: 
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level, format='[%(levelname)s] (%(threadName)-10s) %(message)s')
    threads = []     
    logging.info('Starting JournalBot')
    if args.output != None:
        logging.info('Output: %s' %args.output)
    # Load journal black list: 
    blacklistfile = os.path.join(args.workdir,'blacklist.txt')
    if os.path.exists(blacklistfile):
        f = open(blacklistfile,'r')
        for line in f.readlines():
            blacklist[line.strip()] = True 
    # Load Journal IF db:
    journalfile = os.path.join(args.workdir,'ifactor.txt')
    if os.path.exists(journalfile):
        regex = re.compile(r'^\d+\s([\[A-Za-z&\-\s]+)[ \t]+\d+,?\d+\W+(\d+).\d+')
        with open(journalfile,'r') as f:
            for line in f.readlines()[1:]:
                try: 
                    # 129 Annual Review of Medicine 5,612 12.928 0.0135
                    match = regex.search(line)
                    ifactor[str(match.group(1)).lower()] = int(match.group(2))
                except Exception:
                    print('Could not read line:  %s' % line.strip())

    # Init Posting thread 
    post = threading.Thread(target = post_thread)
    threads.append(post)
    for t in threads: 
        t.start()
    for t in threads: 
        t.join()
        
def _cleanhtml(raw_html):
    cleanr = re.compile('<.*?>')
    cleantext = re.sub(cleanr, '', raw_html.decode('utf-8'))
    cleantext = re.sub(re.compile('&amp;'), '' , cleantext)
    return cleantext

def _statushash(text):
    text = re.sub(re.compile('&amp;'), '' , text)
    text = text.replace('amp', '')
    text = re.sub(re.compile('[^a-zA-Z0-9]'), '' , text)
    return text.lower()

def get_next_paper(paperlist, access_token, CUTOFF, googlesheet=None, googlecred=None):
    Entrez.email = os.environ.get('entrez_email', 'nabil@happykhan.com')
    # Retrieve all papers for each author from file 
    # Use entrez pubmed search
    logging.info('Fetching next paper to post')
    # Retrieve file saved last update time if available 
    search_terms = {} 
    search_file = os.path.join(args.workdir, 'searchterms.txt')
    total_string = ''
    # Handle google sheet as well.
    term_list = [] 
    if googlecred and googlesheet:
        gc = gspread.service_account(filename=googlecred)
        sh = gc.open(googlesheet).sheet1
        term_list = sh.col_values(1)[1:]
    elif os.path.exists(search_file):
        for term in csv.DictReader(open(search_file), delimiter=';'):
            term_list.append(term["search_term"].strip())
    for term in term_list:
        total_string += f' OR ({term})'
        position = round(len(total_string) / 2500)
        if search_terms.get(position):
            search_terms[position] += f" OR ({term})"
        else:
            search_terms[position] = f"({term})"
    entries = [] 
    for search_string in search_terms.values(): 
        try:
            logging.info('Fetching papers relating to ' + search_string)
            handle = Entrez.esearch(db="pubmed", term = search_string, reldate =  CUTOFF) 
            record = Entrez.read(handle)
            handle.close()
            time.sleep(3) 
            if len(record['IdList']) > 0:
                entries += record['IdList']
        except Exception: 
            traceback.print_exc()
            logging.error('Failed to fetch ' + search_string + '. Skipping.')                
    chunks = [entries[x:x+20] for x in range(0, len(entries), 20)]
    for entrylist in chunks:
        result = [] 
        try:
            handle = Entrez.esummary(db="pubmed", id=",".join(entrylist))
            result = Entrez.read(handle)
            handle.close()
            for r in result:
                try:
                    if not blacklist.get(r['FullJournalName'], None): 
                        ifact = ifactor.get(r['FullJournalName'].lower(), 0)
                        newpaper = dict( authors = r['AuthorList'], \
                                            pmid = r['ArticleIds']['pubmed'][0].encode('ascii', 'ignore').decode('utf-8'), \
                                            full_title = _cleanhtml(r['Title'].encode('ascii', 'ignore')), \
                                            journal = r['FullJournalName'], \
                                            ifactor = ifact, \
                                            date = r['EPubDate'])
                        twit_length = 0
                        newpaper['tweet_title'] = newpaper['full_title'][0:(280 - (twit_length + 20))]
                        newpaper['score'] = _calculateScore(newpaper, CUTOFF)
                        newpaper['status_hash'] = _statushash(newpaper['tweet_title'])
                        # Is the paper already there? 
                        if newpaper['score']  > 0 \
                            and newpaper.get('pmid') \
                            and not newpaper['status_hash'] in ( item['status_hash'] for item in paperlist ): 
                            paperlist.append(newpaper)
                except:
                    logging.error('Error in loading %s' %r['Title'].encode('ascii', 'ignore'))
        except:
            logging.error('Error in loading %s' %entrylist[0])                    
    paperlist.sort(key=lambda paper: paper['score'], reverse=True)
    for next_paper in paperlist:
        if next_paper['score'] > 0:
            logging.info(_generate_message(next_paper, access_token, test=True))
    return paperlist[0]

def _generate_message(next_paper, access_token, test=False):
    message = next_paper['tweet_title']
    if test:
        message +=  ' http://test/%s' %next_paper['pmid']
    else:
        message += ' ' + getShortUrl('http://www.ncbi.nlm.nih.gov/pubmed/%s' %next_paper['pmid'], access_token)
    return message

def _calculateScore(paper, CUTOFF):
    dateCutoff = datetime.datetime.combine(datetime.date.today() 
                                           - datetime.timedelta(days=CUTOFF),
                                           datetime.datetime.min.time())
    paper['score'] = 0
    if paper.get('date'):
        try:
            if len(paper['date']) == 8:
                pubDate = datetime.datetime.strptime(paper['date'], "%Y %b")
            else:
                pubDate = datetime.datetime.strptime(paper['date'], "%Y %b %d")
            if pubDate > dateCutoff: 
                paper['score'] = (CUTOFF - (datetime.datetime.today() - pubDate).days) 
                paper['score'] += paper['ifactor'] * 2
        except ValueError:
            return 0
    return paper['score']

def post_thread():
    logging.info('Running posting thread')
    consumer_key = os.environ.get('consumer_key', None)
    consumer_secret = os.environ.get('consumer_secret', None)
    access_token = os.environ.get('access_token', None)
    access_token_secret = os.environ.get('access_token_secret', None)
    bitly_access_token = os.environ.get('bitly_access_token', None)

    if os.path.exists('credentials.json'):
        json_values = json.load(open('credentials.json'))
        consumer_key = json_values['consumer_key']
        consumer_secret = json_values['consumer_secret']
        access_token = json_values['access_token']
        access_token_secret = json_values['access_token_secret']
        bitly_access_token = json_values['bitly_access_token']
        google_sheet_name = json_values.get("google_sheet_name", None)
        google_sheet_credentials = json_values.get("google_sheet_credentials", None)

    if not (consumer_secret and consumer_key and access_token and access_token_secret):
        
        logging.error('TOKENS NOT FOUND - Populate the credentials file, then: source ./credentials')
        sys.exit(-1)

    # OAuth process, using the keys and tokens  
    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)  
    auth.set_access_token(access_token, access_token_secret)  
    
    # Creation of the actual interface, using authentication  
    api = tweepy.API(auth)  
    lastTweetTime = datetime.datetime.utcfromtimestamp(0)
    for status in tweepy.Cursor(api.user_timeline).items(1):
        lastTweetTime = status.created_at 

    #Retrieve existing statuses
    while (True):
        logging.info('Waiting. Last Tweet at ' + str(lastTweetTime))
        if datetime.datetime.now() > (lastTweetTime + datetime.timedelta(minutes=args.tinterval)):
            logging.info('Fetching previous tweets')
            page_list = []
            for page in tweepy.Cursor(api.user_timeline, count=200).pages(10):
                page_list.append(page)
            # Retrieve previous tweets
            title_match = re.compile('(NEW: |OF NOTE: )?(.+) http.+')
            paperlist = []
            for page in (page_list):
                for status in (page):
                    # Tweepy seems to truncate tweets by default. 
                    if status.truncated:
                        clean_status = api.get_status(status.id, tweet_mode='extended')\
                            ._json['full_text'].encode('ascii', 'ignore')
                    else:
                        clean_status = status.text.encode('ascii', 'ignore')
                    clean_status = clean_status.decode('utf-8')

                    if title_match.match(clean_status):
                        tweet_title = title_match.match(clean_status).group(2)
                        status_hash = _statushash(tweet_title)
                        if not status_hash in ( item['status_hash'] for item in paperlist ):
                            paperlist.append(dict(tweet_title = tweet_title, score = 0, posted = True, status_hash=status_hash))
            logging.info('Updated previous tweets')
            logging.info('Last Tweet at ' + str(status.created_at))
            next_paper = get_next_paper(paperlist, bitly_access_token, args.date_cutoff, googlesheet=google_sheet_name, googlecred=google_sheet_credentials)
            if next_paper['score'] > 0:
                try:
                    message = _generate_message(next_paper, bitly_access_token)
                    api.update_status(message)
                    lastTweetTime = datetime.datetime.now()
                except: 
                    logging.error('Error posting status: %s' %message)
            else:
                # All possible posts have score 0, Nothing new to post!
                lastTweetTime = datetime.datetime.now()
        time.sleep(600)

def getShortUrl(longurl, access_token):
    query_params = {'access_token': access_token,
                    'longUrl': longurl} 
    endpoint = 'https://api-ssl.bitly.com/v3/shorten'
    response = requests.get(endpoint, params=query_params, verify=False)
    data = json.loads(response.content)
    return data['data']['url']


if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[1].strip()
        parser = argparse.ArgumentParser(description=desc,epilog=epi)
        parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_argument('--version', action='version', version='%(prog)s ' + meta.__version__)
        parser.add_argument('-o','--output',action='store',help='output prefix')
        parser.add_argument('-t', '--tinterval', action='store',help='Time interval between tweets (minutes)', type=int,default=300)
        parser.add_argument ('--workdir', action='store', help='Working directory', default='config')
        parser.add_argument('-d', '--date_cutoff', action='store',help='Post papers in the last X days [def: 60]', type=int,default=60)	
        args = parser.parse_args()
        if args.verbose: 
            print("Executing @ " + time.asctime())
        main()
        if args.verbose: 
            print("Ended @ " + time.asctime())
        if args.verbose: 
            print('total time in minutes: ' + time.time() - start_time / 60.0)
        sys.exit(0)
    except Exception:
        print('ERROR, UNEXPECTED EXCEPTION')
        traceback.print_exc()
        os._exit(1)       
