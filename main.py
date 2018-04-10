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
        regex = re.compile('^\d+\s([\[A-Za-z&\-\s]+)[ \t]+\d+,?\d+\W+(\d+).\d+')
        with open(journalfile,'r') as f:
            for line in f.readlines()[1:]:
                try: 
                    # 129 Annual Review of Medicine 5,612 12.928 0.0135
                    match = regex.search(line)
                    ifactor[str(match.group(1)).lower()] = int(match.group(2))
                except Exception as e:
                    print 'Could not read line:  %s' % line.strip()

    # Init Posting thread 
    post = threading.Thread(target = post_thread)
    threads.append(post)
    for t in threads: 
        t.start()
    for t in threads: 
        t.join()
        
def _cleanhtml(raw_html):
    cleanr = re.compile('<.*?>')
    cleantext = re.sub(cleanr, '', raw_html)
    return cleantext

def get_next_paper(paperlist, CUTOFF):
    Entrez.email = os.environ.get('entrez_email', 'nabil@happykhan.com')
    # Retrieve all papers for each author from file 
    # Use entrez pubmed search
    logging.info('Fetching next paper to post')
    # Retrieve file saved last update time if available 
    search_terms = {} 
    search_file = os.path.join(args.workdir, 'searchterms.txt')
    if os.path.exists(search_file):
        for term in csv.DictReader(open(search_file), delimiter=';'):
            search_terms[term['search_term']] = term.get('twitter_handle', None)
    for search_term, twit_handle in search_terms.iteritems():
        try:
            logging.info('Fetching papers relating to ' + search_term)
            handle = Entrez.esearch(db="pubmed", term = search_term, reldate =  CUTOFF) 
            record = Entrez.read(handle)
            handle.close()
            time.sleep(3) 
            if len(record['IdList']) > 0:
                entries = record['IdList']
                chunks = [entries[x:x+20] for x in xrange(0, len(entries), 20)]
                for entrylist in chunks:
                    result = [] 
                    try:
                        handle = Entrez.esummary(db="pubmed", id=",".join(entrylist))
                        result = Entrez.read(handle)
                        handle.close()
                    except:
                        logging.error('Error in loading %s' %entrylist[0])
                    for r in result:
                        if not blacklist.get(r['FullJournalName'], None): 
                            ifact = ifactor.get(r['FullJournalName'].lower(), 0)
                            newpaper = dict( authors = r['AuthorList'], \
                                             pmid = str(r['ArticleIds']['pubmed'][0].encode('ascii', 'ignore')), \
                                             full_title = str(_cleanhtml(r['Title'].encode('ascii', 'ignore'))), \
                                             journal = r['FullJournalName'], \
                                             ifactor = ifact, \
                                             date = r['EPubDate'],
                                             twitter_handle = twit_handle)
                            twit_length = 0
                            if newpaper.get('twitter_handle'):
                                twit_length = len(newpaper.get('twitter_handle', ''))
                            newpaper['tweet_title'] = newpaper['full_title'][0:(280 - (twit_length + 30))]
                            newpaper['score'] = _calculateScore(newpaper, CUTOFF)
                            # Is the paper already there? 
                            if newpaper['score']  > 0 \
                               and newpaper.get('pmid') \
                               and not newpaper['tweet_title'] in ( item['tweet_title'] for item in paperlist ): 
                                paperlist.append(newpaper)
        except Exception: 
            traceback.print_exc()
            logging.error('Failed to fetch ' + search_term + '. Skipping.')
    paperlist.sort(key=lambda paper: paper['score'], reverse=True)
    for next_paper in paperlist:
        if next_paper['score'] > 0:
            logging.info(_generate_message(next_paper, test=True))
    return paperlist[0]

def _generate_message(next_paper, test=False):
    message = next_paper['tweet_title']
    if test:
        message +=  ' http://test/%s' %next_paper['pmid']
    else:
        message += ' ' + getShortUrl('http://www.ncbi.nlm.nih.gov/pubmed/%s' %next_paper['pmid'])
    if next_paper.get('twitter_handle'):
        join_word = ''
        if next_paper.get('twitter_handle')[0] == '@':
            join_word = 'by'
        message += ' by %s' %next_paper.get('twitter_handle')
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
            page_list = []
            for page in tweepy.Cursor(api.user_timeline, count=200).pages(10):
                page_list.append(page)
            # Retrieve previous tweets
            title_match = re.compile('(NEW: |OF NOTE: )?(.+) http.+')
            paperlist = []
            for page in reversed(page_list):
                for status in reversed(page):
                    clean_status = status.text.encode('ascii', 'ignore')
                    if title_match.match(clean_status):
                        tweet_title = title_match.match(clean_status).group(2)
                        if not tweet_title in ( item['tweet_title'] for item in paperlist ):
                            paperlist.append(dict(tweet_title = tweet_title, score = 0, posted = True))
            logging.info('Last Tweet at ' + str(status.created_at))
            next_paper = get_next_paper(paperlist, args.date_cutoff)
            if next_paper['score'] > 0:
                try:
                    api.update_status(_generate_message(next_paper))
                    lastTweetTime = datetime.datetime.now()
                except: 
                    logging.error('Error posting status: %s' %message)
            else:
                # All possible posts have score 0, Nothing new to post!
                lastTweetTime = datetime.datetime.now()
        time.sleep(600)

def getShortUrl(longurl):
    query_params = {'access_token': '057b8e5650460b576694ead637996aec16c405f3',
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
        parser.add_argument ('workdir', action='store', help='Working directory')
        parser.add_argument('-d', '--date_cutoff', action='store',help='Post papers in the last X days [def: 30]', type=int,default=30)	
        args = parser.parse_args()
        if args.verbose: 
            print "Executing @ " + time.asctime()
        main()
        if args.verbose: 
            print "Ended @ " + time.asctime()
        if args.verbose: 
            print 'total time in minutes:',
        if args.verbose: 
            print (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)       
