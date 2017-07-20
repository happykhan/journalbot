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

TODO: 
    * Post with author as twitter handle.
"""


class Paper(object):
    def __init__(self, citationCount, authors, mainAuthor, pubmedid, title,
                 date,  journal, ifactor=2, score=0, posted= False):
        self.citationCount = citationCount
        self.authors = authors
        self.mainAuthor = mainAuthor
        self.pubmedid = pubmedid
        self.title = title
        self.date = date
        self.posted = posted
        self.score = score
        self.journal = journal
        self.ifactor = ifactor
        

    def isOfNote(self):
        CUTOFF = 20
        import datetime        
        pubDate = datetime.datetime.strptime(self.date, "%Y/%m/%d %H:%M")                    
        yearcount = ( datetime.date.today().year - pubDate.year) 
        if yearcount < 1: yearcount = 1
        citsperyear = int(round(int(self.citationCount) / yearcount,0))     
        if citsperyear > CUTOFF: return True
        else: return False

    def isNew(self):
        import datetime
        CUTOFF = 10
        pubDate = datetime.datetime.strptime(self.date, "%Y/%m/%d %H:%M")          
        dateCutoff = datetime.datetime.combine(datetime.date.today() - datetime.timedelta(days=CUTOFF), datetime.datetime.min.time())     
        if pubDate > dateCutoff: return True
        else: return False

    def calculateScore(self):
        import datetime
        import random
        CUTOFF = 180
        dateCutoff = datetime.datetime.combine(datetime.date.today() - datetime.timedelta(days=CUTOFF), datetime.datetime.min.time())
        pubDate = datetime.datetime.strptime(self.date, "%Y/%m/%d %H:%M")
        if pubDate > dateCutoff: 
            self.score = (CUTOFF - (datetime.datetime.today() - pubDate).days) 
        else: 
            self.score = 0
        yearcount = ( datetime.date.today().year - pubDate.year) -1 
        if yearcount < 1: 
            yearcount = 1
        citsperyear = int(round(int(self.citationCount) / yearcount,0))
        if len(self.authors) > 5: 
            if self.mainAuthor in self.authors[2:-2]:
                self.score = self.score -  (len(self.authors)  * 5 )
#        self.score += citsperyear
        self.score += (self.ifactor) 
        self.score += random.randint(1,5)

    def __str__(self):
        return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%d' %(self.mainAuthor,  self.title, self.authors, self.date, self.citationCount, self.pubmedid, self.posted, self.score, self.journal, self.ifactor)

    @classmethod
    def from_string(cls, string): 
        linArray = string.strip().split('\t') 
    
        return cls(linArray[4], ast.literal_eval(linArray[2]), linArray[0], linArray[5], linArray[1], linArray[3], linArray[8],score=int(linArray[7]),ifactor =int(linArray[9]),posted=linArray[6]   )


import sys, os, traceback, argparse, warnings
import time, re
import __init__ as meta 
import threading
import logging
from time import gmtime
import ast 

epi = "Licence: "+meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"

paperlist = []
lock = threading.Lock()
ifactor = {} 
blacklist = {} 

consumer_key = None
consumer_secret = None
access_token =  None
access_token_secret = None

def main ():

    global args
    print 'Starting JournalBot'
    if args.output != None:
        print 'Output: ' + args.output
    global consumer_key
    global consumer_secret
    global access_token
    global access_token_secret
    consumer_key = os.environ.get('consumer_key', None)
    consumer_secret = os.environ.get('consumer_secret', None)
    access_token = os.environ.get('access_token', None)
    access_token_secret = os.environ.get('access_token_secret', None)
    if not (consumer_secret and consumer_key and access_token and access_token_secret):
        logging.error('TOKENS NOT FOUND - Populate the credentials file, then: source ./credentials')
        sys.exit(0)
    # Load journal black list: 
    blacklistfile = os.path.join(args.workdir,'blacklist.txt')
    if os.path.exists(blacklistfile):
        f = open(blacklistfile,'r')
        for line in f.readlines():
            blacklist[line.strip()] = True 
    
    # Load up Papers db. 
    paperlistfile = os.path.join(args.workdir,'paperlist.txt')
    if os.path.exists(paperlistfile):
        f = open(paperlistfile,'r')
        for line in f.readlines():
            newpaper = Paper.from_string(line)
            newpaper.calculateScore()
            paperlist.append(newpaper)     
    # Load Journal db:
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
    

    # Init paper db update thread 
    if args.verbose: logging.basicConfig(level=logging.DEBUG,format='[%(levelname)s] (%(threadName)-10s) %(message)s',)
    else:     logging.basicConfig(level=logging.INFO,format='[%(levelname)s] (%(threadName)-10s) %(message)s',)
    threads = [] 
    # Init citation count thread
    # Init Posting thread 
    update = threading.Thread(target=updateThread,args=(args.updatehours,))
    threads.append(update)
    post = threading.Thread(target=postThread, args=(paperlist,))
    threads.append(post)
    for t in threads: t.start()
    for t in threads: t.join()



def updateThread(updatehours):
    from Bio import Entrez
    Entrez.email = os.environ.get('entrez_email', None)
    # Retrieve all papers for each author from file 
    # Use entrez pubmed search
    global lock
    logging.info('Launching update thread')
    updatefile = os.path.join(args.workdir,'lastupdate.txt')
    paperlistfile = os.path.join(args.workdir,'paperlist.txt')
    # Retrieve file saved last update time if available 
    while (True):
        lastupdate = 0        
        if os.path.exists(updatefile):    
            try:
                lastupdate = float(open(updatefile).readline().strip())
            except:
                logging.error('Can not get last update time. File corrupted')        
        logging.info('Last updated %s' %time.ctime(lastupdate))
        if lastupdate + (updatehours * 3600) < time.time():
            logging.info('Running update cycle')
            lock.acquire()
            global paperlist            
            authorlist = [] 
            f = open(os.path.join(args.workdir,'authorlist.txt'))
            for line in f.readlines():
                authorlist.append(line.strip())                
            for mainAuthor in authorlist:
                count = 0 
                while count < 3: 
                    try:
                        logging.info('Fetching papers from ' + mainAuthor)                        
                        handle = Entrez.esearch(db="pubmed", term=mainAuthor, field='author')
                        record = Entrez.read(handle)
                        time.sleep(3) 
                        handle.close()
                        entries = record['IdList']
			print 'Retrieved the following entries:'
			for entry in entries : print entry
                        entrylist = []
                        count = 0
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
                                    ifact = ifactor.get(r['FullJournalName'].lower(), 2)
                                    newpaper = Paper(r['PmcRefCount'], str(r['AuthorList']).encode('ascii', 'ignore'), \
                                                         mainAuthor, str(r['ArticleIds']['pubmed'][0].encode('ascii', 'ignore')), \
                                                         str(r['Title'].encode('ascii', 'ignore')), \
                                                         str(r['History']['pubmed'][0].encode('ascii', 'ignore')), r['FullJournalName'], ifact) 
                                    newpaper.calculateScore()
                                    # Is the paper already there? 
                                    duplicate = False
                                    for oldpaper in paperlist:
                                        if newpaper.title == oldpaper.title:
                                            duplicate = True
                                            oldpaper.citationCount = newpaper.citationCount
                                            oldpaper.calculateScore()
                                    if not duplicate: 
                                        paperlist.append(newpaper)
                        break
                    except Exception: 
                        count += 1 
                        logging.error('Failed to call authors update' + mainAuthor + '. Waiting 30 seconds.')
                        time.sleep(30) 
                        #continue
            paperlist.sort(key=lambda paper: (paper.score, paper.citationCount), reverse=True)
            f = open(paperlistfile, 'w')            
            for paper in paperlist:
                f.write('%s\n' %paper)
            f.close()
            lock.release()
            logging.info('Paper database now with %d records' %len(paperlist))            
            # update sucessful run          
            lastupdate = time.time()
            logging.info('Finished update cycle at %s' %time.ctime(lastupdate))
            f = open(updatefile, 'w')
            f.write(str(lastupdate)+'\n')
            f.close()
        time.sleep(600)      

def exportPaperDb():
    print 'not done'


def postThread(initPaperlist):
    time.sleep(5)              
    logging.info('Running posting thread')
    import tweepy  
    import random       
    import datetime          
    global consumer_key
    global consumer_secret
    global access_token
    global access_token_secret
    # OAuth process, using the keys and tokens  
    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)  
    auth.set_access_token(access_token, access_token_secret)  
    # Creation of the actual interface, using authentication  
    api = tweepy.API(auth)  
    paperlistfile = os.path.join(args.workdir,'paperlist.txt')
    lastTweetTime = datetime.datetime.utcfromtimestamp(0)
    for status in tweepy.Cursor(api.user_timeline).items(1):
        lastTweetTime = status.created_at 
    #Retrieve existing statuses
    page_list = [] 
    for page in tweepy.Cursor(api.user_timeline, count=200).pages(16):
        page_list.append(page)
    # Intial check of paper post status.
    for paper in initPaperlist:
        message = ''
        if paper.isNew() == True: message += 'NEW: '
        if paper.isOfNote() == True: message += 'OF NOTE: '
        tweetTitle = paper.title[0:(140-(len(message)+26))]          
        for page in reversed(page_list):
            for status in reversed(page):
                if status.text.encode('ascii', 'ignore').find(tweetTitle) != -1: 
                    paper.posted = 'True'       
    f = open(paperlistfile, 'w')            
    initPaperlist.sort(key=lambda paper: (paper.score, paper.citationCount), reverse=True)
    for paper in initPaperlist:
        if not blacklist.get(paper.journal, None):         
            f.write('%s\n' %paper)
    f.close()
    while (True):
        if not lock.locked():   
            logging.info('Last Tweet at %s' %status.created_at)                
            if datetime.datetime.now() > (lastTweetTime + datetime.timedelta(minutes=300)):
                logging.debug('Posting: %s' %status.created_at)
                global paperlist                
                for paper in paperlist:
                    message = ''
                    if paper.isNew() == True: message += 'NEW: '
                    if paper.isOfNote() == True: message += 'OF NOTE: '
                    tweetTitle = paper.title[0:(140-(len(message)+26))]                              
                    # Check if posted before from tweetlist
                    for page in reversed(page_list):
                        for status in reversed(page):
                            if status.text.encode('ascii', 'ignore').find(tweetTitle) != -1: 
                                paper.posted = 'True'                              
                    if paper.posted == 'False':
                        message += tweetTitle
                        message += ' ' + getShortUrl('http://www.ncbi.nlm.nih.gov/pubmed/%s' %paper.pubmedid)
                        try:                      
                            api.update_status(message)
                            lastTweetTime = datetime.datetime.now()
                            paper.posted = 'True'                                      
                        except: 
                            logging.error('Error posting status: %s' %message)
                        f = open(paperlistfile, 'w')            
                        for paper in paperlist:
                            f.write('%s\n' %paper)
                        f.close()             
                        break; 
        time.sleep(600)               

def getShortUrl(longurl):            
    import requests
    import json

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
        parser.add_argument('-u','--updatehours',action='store',help='interval of hours to run update',type=int,default=4)        
        parser.add_argument ('workdir', action='store', help='Working directory')
        args = parser.parse_args()
        if args.verbose: print "Executing @ " + time.asctime()
        main()
        if args.verbose: print "Ended @ " + time.asctime()
        if args.verbose: print 'total time in minutes:',
        if args.verbose: print (time.time() - start_time) / 60.0
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
