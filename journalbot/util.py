import re
import datetime
import os

def cleanhtml(raw_html):
    cleanr = re.compile("<.*?>")
    cleantext = re.sub(cleanr, "", raw_html.decode("utf-8"))
    cleantext = re.sub(re.compile("&amp;"), "", cleantext)
    return cleantext


def statushash(text):
    text = re.sub(re.compile("&amp;"), "", text)
    text = text.replace("amp", "")
    text = re.sub(re.compile("[^a-zA-Z0-9]"), "", text)
    return text.lower()


def generate_message(next_paper):
    message = next_paper["tweet_title"]
    message += " " + "http://www.ncbi.nlm.nih.gov/pubmed/%s" % next_paper["pmid"]
    return message


def load_credentials():
    consumer_key = os.environ.get("consumer_key", None)
    consumer_secret = os.environ.get("consumer_secret", None)
    access_token = os.environ.get("access_token", None)
    access_token_secret = os.environ.get("access_token_secret", None)
    return consumer_key, consumer_secret, access_token, access_token_secret


def calculateScore(paper, CUTOFF):
    dateCutoff = datetime.datetime.combine(
        datetime.date.today() - datetime.timedelta(days=CUTOFF),
        datetime.datetime.min.time(),
    )
    boost_date = datetime.datetime.combine(
        datetime.date.today() - datetime.timedelta(days=round(CUTOFF/4)),
        datetime.datetime.min.time(),
    )
    paper["score"] = 0
    if paper.get("date"):
        try:
            if len(paper["date"]) == 8:
                pubDate = datetime.datetime.strptime(paper["date"], "%Y %b")
            else:
                pubDate = datetime.datetime.strptime(paper["date"], "%Y %b %d")
            if pubDate > dateCutoff:
                paper["score"] = 100
                paper["score"] += (pubDate - dateCutoff).days * 2 
                if pubDate > boost_date :
                    paper["score"] += 200
                if int(paper["pub_count"]) > 50:
                    paper["score"] += 100
                if int(paper["pub_count"]) < 10:
                    paper["score"] = paper["score"] - 150
        except ValueError:
            return 0
    return paper["score"]
