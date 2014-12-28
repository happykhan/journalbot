def citationThread():
    # Retrieve citation count using scholar.py    
    logging.debug('running citation')
    global paperlist
    while (True):
          time.sleep(1)      
          logging.debug('citation alive')
          import subprocess
          for curpaper in paperlist: 
               if curpaper.citationCount == '-1':
                    import random
                    time.sleep(random.randint(10,90))
                    print 'Z%sZ' %curpaper.mainAuthor
                    print 'Z%sZ'%curpaper.title
                    process = subprocess.Popen([sys.executable, "scholar.py",'--csv', curpaper.mainAuthor, '--title-only', '--all', curpaper.title ], stdout=subprocess.PIPE,shell = True)
                    result = process.communicate()[0]
                    print 'zz %s' %result
  