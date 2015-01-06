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
"""


class Paper(object):
    def __init__(self, citationCount, authors, mainAuthor, pubmedid, title,
                 date, posted, score):
        self.citationCount = citationCount
        self.authors = authors
        self.mainAuthor = mainAuthor
        self.pubmedid = pubmedid
        self.title = title
        self.date = date
        self.posted = posted
        self.score = score

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
        CUTOFF = 28
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
        if pubDate > dateCutoff: self.score = (CUTOFF - (datetime.datetime.today() - pubDate).days) 
        else: self.score = 0
        yearcount = ( datetime.date.today().year - pubDate.year) -1 
        if yearcount < 1: yearcount = 1
        citsperyear = int(round(int(self.citationCount) / yearcount,0))
        self.score += citsperyear
        self.score += random.randint(1,10)

    def __str__(self):
        return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d' %(self.mainAuthor,  self.title, self.authors, self.date, self.citationCount, self.pubmedid, self.posted,self.score)



import sys, os, traceback, argparse
import time
import __init__ as meta 
import threading
import logging
from time import gmtime

epi = "Licence: "+meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"

paperlist = []
lock = threading.Lock()

def main ():

    global args
    print 'Hello world!'
    if args.output != None:
        print 'Output: ' + args.output

    # Load up Papers db. 
    paperlistfile = os.path.join(args.workdir,'paperlist.txt')
    if os.path.exists(paperlistfile):
        f = open(paperlistfile,'r')
        for line in f.readlines():
            linArray = line.strip().split('\t') 
            newpaper = Paper(linArray[4], linArray[2], linArray[0], linArray[5], 
                             linArray[1], linArray[3], linArray[6],int(linArray[7]))
            paperlist.append(newpaper)     

    # Init paper db update thread 
    if args.verbose: logging.basicConfig(level=logging.INFO,format='[%(levelname)s] (%(threadName)-10s) %(message)s',)
    else:     logging.basicConfig(level=logging.DEBUG,format='[%(levelname)s] (%(threadName)-10s) %(message)s',)
    threads = [] 
    # Init citation count thread
    # Init Posting thread 
    update = threading.Thread(target=updateThread,args=(args.updatehours,))
    threads.append(update)
    post = threading.Thread(target=postThread)
    threads.append(post)
    for t in threads: t.start()
    for t in threads: t.join()



def updateThread(updatehours):
    from Bio import Entrez
    Entrez.email = "nabil@happykhan.com"
    # Retrieve all papers for each author from file 
    # Use entrez pubmed search
    global paperlist, lock
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
            authorlist = [] 
            f = open(os.path.join(args.workdir,'authorlist.txt'))
            for line in f.readlines():
                authorlist.append(line.strip())                
            for mainAuthor in authorlist:
                handle = Entrez.esearch(db="pubmed", term="%s [AU]" %mainAuthor,retmax=10000)
                record = Entrez.read(handle)
                handle.close()
                entries = record['IdList']
                entrylist = []
                count = 0
                for e in entries:
                    entrylist.append(e)
                    count += 1
                    if count % 20 == 0 or count == len(entries): 
                        result = [] 
                        try:
                            handle = Entrez.esummary(db="pubmed", id=",".join(entrylist))
                            result = Entrez.read(handle)
                            handle.close()
                        except:
                            logging.error('Error in loading %s' %entrylist[0])
                        entrylist = []                                      
                        for r in result:
                            newpaper = Paper(r['PmcRefCount'], str(r['AuthorList']).encode('ascii', 'ignore'), \
                                             mainAuthor, str(r['ArticleIds']['pubmed'][0].encode('ascii', 'ignore')), \
                                             str(r['Title'].encode('ascii', 'ignore')), \
                                             str(r['History']['pubmed'][0].encode('ascii', 'ignore')), False,0) 
                            newpaper.calculateScore()
                            # Is the paper already there? 
                            duplicate = False
                            for oldpaper in paperlist:
                                if newpaper.title == oldpaper.title:
                                    duplicate = True
                                    oldpaper.citationCount = newpaper.citationCount
                                    oldpaper.calculateScore()
                            if not duplicate: paperlist.append(newpaper)
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


def postThread():
    logging.info('Running posting thread')
    global paperlist
    import tweepy  
    import random       
    import datetime          
    # Consumer keys and access tokens, used for OAuth  
    consumer_key = 'IWz0nMsnkbX1hNnxxQMByCdIU'  
    consumer_secret = 'lvgfYnyct8GorAbRrkCuznRFbNOTSWpghiTLHAdanPTJuvwGCs'  
    access_token = '2938273641-7Uo4d2TpP8fPWtFGRkNJsnUqH1ZauQq5n9ugVtq'  
    access_token_secret = '8yrFOztLdrmPIhtv6jcBFzA7Z9fNqNAnQ3szIDFtQ8Hdi'  

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
    for paper in paperlist:
        message = ''
        if paper.isNew() == True: message += 'NEW: '
        if paper.isOfNote() == True: message += 'OF NOTE: '
        tweetTitle = paper.title[0:(140-(len(message)+26))]          
        for page in reversed(page_list):
            for status in reversed(page):
                if status.text.encode('ascii', 'ignore').find(tweetTitle) != -1: 
                    paper.posted = 'True'       
    f = open(paperlistfile, 'w')            
    for paper in paperlist:
        f.write('%s\n' %paper)
    f.close()                 
    while (True):
        if not lock.locked():   
            logging.info('Last Tweet at %s' %status.created_at)                
            if datetime.datetime.now() > (lastTweetTime + datetime.timedelta(minutes=55)):    
                logging.debug('Posting: %s' %status.created_at)
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
        parser.add_argument('-u','--updatehours',action='store',help='interval of hours to run update',type=int,default=3)        
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
