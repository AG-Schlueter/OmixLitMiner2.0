try:
    import google.colab
    # Running on Google Colab, so install Biopython first
    !pip install biopython
except ImportError:
    pass
import re
import urllib
import numpy as np
import matplotlib.pyplot as plt
import math
import time
from tqdm import tqdm
import datetime
from datetime import date
from collections import OrderedDict, defaultdict
from nltk.probability import FreqDist
from nltk.tokenize import word_tokenize
import Bio
from Bio import Entrez, Medline
import pandas as pd
import plotly.express as px
import io
import requests
from openpyxl import Workbook
import ipywidgets as widgets
from google.colab import files
import nltk
#nltk.download('punkt')

import copy
#@title Install dependencies and enter input

Mail = 'an.gocke@uke.de'

#@markdown  Specify whether you want to search PubMed directly or whether you want to use a more generalised approach, utilising the AI-based PubTator3.0. For a list of the differences, please check the publication ().
PubMed = True #@param {type:"boolean"}
PubTator3 = False #@param {type:"boolean"}

#@markdown  Specify the title of the output file. If you have searched in the title only, it will be added to the filename
Filename = 'Filename'#@param {type:"string"}

Protein_List = 'ALB' #@param {type:"string"}
#@markdown  Specify the type of the Input for the Protein List, set a tick either for the Accession or the Gene Name
Accession = False #@param {type:"boolean"}
Gene_name = True #@param {type:"boolean"}
Gene_name = True

#if (Accession and Gene_name) or (not Accession and not Gene_name):
#  raise RuntimeError(f'You need to enter either Accession or Gene_name')

#@markdown  Taxonomy ID (e.g. homo sapiens: 9606; mus musculus: 10090)
TaxID = "9606"#@param {type:"string"}

#@markdown  If you want to enter multiple keywords, please separate them by semicolons.
Keywords = 'cancer'#@param {type:"string"}

#@markdown  Checks if the keywords are mentioned  in the title only or in the title and the abstract
KeywordInTitleOnly = True #@param {type:"boolean"}
if KeywordInTitleOnly == True:
  Title = True
  Title_Abstract = False
else:
  Title = False
  Title_Abstract = True

#@markdown  Set number of how many papers should be maximally retrieved (max 1000)
maxPaper = 1000 #@param {type:"integer"}

if maxPaper > 1000:
  maxPaper = 1000

#@markdown  Specify the order for the publication output - standard is relevance
Publication_date = False #@param {type:"boolean"}
Relevance = True #@param {type:"boolean"}
if Publication_date and Relevance:
  raise RuntimeError(f'You can only sort by the publication date or the relevance')


if Publication_date:
  order = "pub_date"
else:
  order = "relevance"


protList = Protein_List.split()
keywords = Keywords

if Accession:
  IDType = "Accession"
else:
  IDType = "Gene"

def requestOnUniProtWebpage(uniprotID, idtype, taxid, proteinrequired):
  start = "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true"
  middle = '+AND+accession:'\
          if str(idtype).lower() == 'accession' \
          else '{}{}{}'.format("+AND+organism_id:", taxid, "+AND+gene:")
  end = '&format=tsv'
  contents = urllib.request.urlopen('{}{}{}{}'.format(start,
                 middle, uniprotID, end)).read().decode("utf-8")

  data = io.StringIO(contents)
  df = pd.read_csv(data, sep="\t", header=0)
  return df

def getUniProtSynonyms(uniprotID, idtype, taxid):
  contents = requestOnUniProtWebpage(uniprotID, idtype, taxid, False)
  if len(contents['Gene Names']) > 0:
    syns = [str(x) for x in contents['Gene Names']][0]
    syns = re.sub(r'\s', ',', syns)
    syns = syns.split(",")
    syns.append(uniprotID)
  else:
    syns = [uniprotID]
  if len(contents['Protein names']) > 0:
    prots = ([str(x) for x in contents['Protein names']][0])
    split_parts = re.split(r"\s(?=\()", prots, maxsplit=1)
    main_part = split_parts[0].strip().replace(' ', '-')
    remaining_part = prots[len(main_part):]
    parentheses_parts = []
    stack = []
    current = ""
    for char in remaining_part:
        if char == "(":
            if stack:
                current += char
            stack.append(char)
        elif char == ")":
            stack.pop()
            if not stack:
                parentheses_parts.append(current.strip().replace(' ', '-'))
                current = ""
            else:
                current += char
        elif stack:
            current += char
    contents = [main_part] + parentheses_parts
  else:
    contents = ''
  return syns, contents


def uniProtQuery(uniprotID, idtype, taxid, keyword, tiab, mail, maxdate):

  synonyms, protnames = getUniProtSynonyms(uniprotID, idtype, taxid)
  searchsummary = OrderedDict()
  searchsummary['uniprotID'] = uniprotID
  searchsummary['idtype'] = idtype
  searchsummary['taxid'] = taxid
  searchsummary['synonyms'] = synonyms
  searchsummary['protein_names'] = protnames
  searchsummary['keywords'] = keyword
  searchsummary['totalResults'] = 0
  searchsummary['reviewResults'] = 0
  searchsummary['category'] = 0
  searchsummary['false'] = 0
  uniprotSyn = [uniprotID]

  papercol = Papercollection(uniprotID, idtype, taxid, synonyms, protnames,\
                             mail, KeywordInTitleOnly, keyword,  maxdate)
  searchsummary['totalResults'] = papercol.resultcnt()
  searchsummary['reviewResults'] = papercol.revcnt()
  if searchsummary['reviewResults'] > 0:
    searchsummary['category'] = 1
  elif searchsummary['totalResults'] > 0:
    searchsummary['category'] = 2
  elif searchsummary['synonyms'] == uniprotSyn:
    searchsummary['category'] = 4
    searchsummary['synonyms'] = []
  else:
    searchsummary['category'] = 3
  return searchsummary, papercol, papercol.resultcnt(), papercol.revcnt()


class Papercollection:
  def __init__(self, uniprotID, idtype, taxid, synonyms, protnames,
                mail, titleonly, keyword, maxdate):
    if PubMed == True:
      self._titlesPapercollection = ['Titles','Abstracts','Years','Authors',\
                             'Affiliation','Country','PMID','PTyp','Reviews',\
                              'Category','SearchedFor']
      self._shortTitles = ['TI','AB','DP','AU','AD','PL','PMID','PT']
      self._uniprotID = uniprotID
      self._idtype = idtype
      self._taxid = taxid
      self._titleonly = titleonly
      self._mail = mail
      self._syns = synonyms
      self._prots = protnames
      self._papercollection = OrderedDict()
      self._PMIDdict = dict()
      self._maxdate = maxdate
      self._keywordlist = keyword.split(';')
      self._resultcount = 0
      self._reviewcount = 0
      fields = '[TI]' if titleonly else ''
      if len(self._keywordlist) == 1:
        self._keyword = '{}{}'.format(keyword, fields)

      #mehrere keywords durchgehen
      elif len(self._keywordlist) >1:
        for idx in range(len(self._keywordlist)):
          self._keywordlist[idx] = '{}{}'.format(self._keywordlist[idx],fields)
        if Together == True:
          placeholder = ' AND '.join(self._keywordlist)
        else:
          placeholder = ' OR '.join(self._keywordlist)
        self._keyword = '({})'.format(placeholder)

      for title in self._titlesPapercollection:
        self._papercollection[title] = list()
      records, term = self.searchAll(self._syns, self._prots)
      self.workThroughRecords(records, term, True)

    elif PubTator3 == True:
      self._titlesPapercollection = ['Titles','Authors','Journal','Years',\
                             'PMID','DOI','Link']
      self._syns = synonyms
      self._prots = protnames
      self._papercollection = OrderedDict()
      self._keywordlist = keyword.split(';')
      self._keywordlist = list(self._keywordlist)

      if len(self._keywordlist) == 1:
        self._keyword = self._keywordlist
      #mehrere keywords durchgehen
      elif len(self._keywordlist) >1:
        if Together == True:
          placeholder = '%20AND%20'.join(self._keywordlist)
          self._keyword = '{}'.format(placeholder)
        else:
          self._keyword = self._keywordlist
      for title in self._titlesPapercollection:
        self._papercollection[title] = list()

      records, term = self.requestPubTator(self._syns, self._prots)
      self.parseResultsOPM(records)

  def revcnt(self):
    return(self._reviewcount)

  def categoryOnePaper(self, paper, synonym, keywords, tion):
    allwords, keys = set(), set()
    category = 3
    review = False
    if not {'TI','AB','PT'}.issubset(set(paper.keys())) or len(paper['TI']) <3:
      return review, category
    ti = paper['TI']
    if not tion:
      ti = '{} {}'.format(ti, paper['AB'])
    ti = re.sub(r'\.|\:|\,|\[|\]|\)|\(','', ti)
    ti = re.sub(r'\-',' ', ti).split()
    #set comprehension
    allwords = {word.lower() for word in ti}
    for word in keywords:
      oneSet = word.split(' ')
      for oneterm in oneSet:
        if len(oneterm) > 0:
          keys.add(oneterm.lower())
    ptyp = {onetype.lower() for onetype in paper['PT']}
    if not keys.isdisjoint(allwords):
      category = 2
      reviewtrue = re.search('review', ''.join(ptyp), flags=re.IGNORECASE)
      if reviewtrue:
        category = 1
        review = True
        self._reviewcount += 1
    else:
      category = 3
    return review, category

  def categoryOnePaperMoreSyns(self, paper, synonyms, keywords, tion):
    reBest, catBest = False, 4
    for syn in synonyms:
      re, cat = self.categoryOnePaper(paper, syn, keywords, tion)
      if not reBest:
        reBest = re
      if cat >0 and cat < catBest:
        catBest = cat
    return reBest, catBest

  def workThroughRecords(self, records, syn, multiple):
    for paper in records:
      if 'DP' in paper:
        year_re = re.search(r'\b\d{4}\b', paper['DP'])
        year = int(year_re.group())
        paper['DP'] = year
        self.onePaperRecord(paper)
        self._papercollection['SearchedFor'].append(self._termquery)
      if 'PMID' in paper.keys():
        if paper['PMID'] in self._PMIDdict:
          continue
        else:
          self._PMIDdict[paper['PMID']] = True

      review, category = False, 0
      if len(paper.keys())>2:
        review, category = self.categoryOnePaperMoreSyns(paper, syn,\
                                                       self._keywordlist, self._titleonly)
      self._papercollection['Category'].append(category)
      self._papercollection['Reviews'].append(review)

  def getPapercolletion(self):
    return self._papercollection

  def resultcnt(self):
    if PubMed == True:
      if self._papercollection['Titles'] == ['nan']:
        x = 0
      else:
        x = len(self._papercollection['Titles'])
    else:
      x = self._resultcount
    return x

  def fetch_handle(self, id_list):
    Entrez.email = self._mail
    handle = Entrez.efetch(db='pubmed', id=id_list, rettype='medline',
                                                    retmode='text')
    return(Medline.parse(handle))

  def searchAll(self,synonyms, proteinnames):
    syns = []
    if len(synonyms) > 0 or len(proteinnames) > 0:
      syns.extend(synonyms)
      syns.extend(proteinnames)

    syns.append(self._uniprotID)
    time.sleep(0.3)
    if len(syns) > 0:
      results, term = self.searchAllSyns(syns)
      id_list = results['IdList']
    return(self.fetch_handle(id_list), term)

  def searchAllSyns(self, syns):
    syns = [syn.replace("(", "") for syn in syns]
    syns = [syn.replace(")", "") for syn in syns]
    syns = [syn.replace(",", "%2C") for syn in syns]
    self._termquery ='{} AND ({}[TI])'.format(self._keyword, syns[0])
    for syn in syns[1:-1]:
      self._termquery = re.sub(r'\)', '',self._termquery)
      self._termquery = '{} OR {}[TI])'.format(self._termquery, syn)
    Entrez.email = self._mail
    handle = Entrez.esearch(db='pubmed',
                            sort=order,
                            retmax=maxPaper,
                            retmode='xml',
                            term=self._termquery,
                            maxdate=self._maxdate,
                            mindate='1800/01/01')
    results = Entrez.read(handle)
    handle.close()
    return results, self._termquery

  def onePaperRecord(self, paper):
    for idx, title in enumerate(self._shortTitles):
      if title in paper:
        self._papercollection[self._titlesPapercollection[idx]].append(paper[title])
      else:
        self._papercollection[self._titlesPapercollection[idx]].append('nan')

  def requestPubTator(self, synonyms, proteinnames):
    start = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text="
    pubtator = "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text="
    combAnd = '%20AND%20'
    syns = []
    if len(synonyms) > 0 or len(proteinnames) > 0:
      if isinstance(synonyms, str):
        synonyms = [synonyms]
      if isinstance(proteinnames, str):
        proteinnames = [proteinnames]
      syns =  synonyms + proteinnames
    synquery = []
    if len(self._keyword) >1:
      for keyword in self._keyword:
          for syn in syns:
            if isinstance(syn, list):
              syn = " ".join(syn)
            syn = re.sub(r"[()/]", "", syn)
            syn = '({}{}{})'.format(keyword, combAnd, syn)
            synquery.append(syn)
    else:
        for syn in syns:
          if isinstance(syn, list):
            syn = " ".join(syn)
          syn = re.sub(r"[()/]", "", syn)
          syn = '({}{}{})'.format(self._keyword, combAnd, syn)
          synquery.append(syn)
    protQuery = '%20OR%20'.join(synquery)
    time.sleep(0.3)

    url = '{}{}'.format(start,  protQuery)
    pubLink = '{}{}'.format(pubtator,  protQuery)
    if Title:
      title = "&sections=title"
      url = '{}{}'.format(url, title)
      pubLink = '{}{}'.format(pubLink, title)
    elif Title_Abstract:
      title_abstract = "&sections=title,abstract"
      url = '{}{}'.format(url, title_abstract)
      pubLink = '{}{}'.format(pubLink, title_abstract)
    if Publication_date:
      pub = "&sort=date%20desc"
      url = '{}{}'.format(url, pub)
      pubLink = '{}{}'.format(pubLink, pub)

    contents = requests.get(url)
    results = contents.json()
    self._url = pubLink
    #print(url)
    return results, url

  def parseResultsOPM(self, records):
    self._resultcount = 0
    self._reviewcount = 0

    results = records["results"]
    for result in results:
      self._papercollection['Titles'].append(result["title"])
      if "authors" in result:
        self._papercollection['Authors'].append(result["authors"])
      else:
        self._papercollection['Authors'].append("NaN")
      self._papercollection['PMID'].append(result["pmid"])
      self._papercollection['Journal'].append(result["journal"])
      if "doi" in result:
        self._papercollection['DOI'].append(result["doi"])
      else:
        self._papercollection['DOI'].append("NaN")
      date = result["meta_date_publication"]
      date_re = re.search(r'\b\d{4}\b', date)
      year = int(date_re.group())
      self._papercollection['Years'].append(year)
      self._papercollection['Link'].append(self._url)

    facets = records["facets"]
    if "facet_fields" in facets:
      facet_fields =facets["facet_fields"]
      facet_types = facet_fields["type"]
      for ftype in facet_types:
        if ftype["name"] == "Review":
          self._reviewcount =ftype["value"]
    else:
      self._reviewcount = 0

    self._resultcount = records["count"]

def runOLM():
  #OmixLitMiner: 13.03.2019
  #maxdate = datetime.date(2019,3,13)
  maxdate= date.today().strftime('%Y/%m/%d')
  papercol_titles = ['Titles','Abstracts','Years','Authors', 'Affiliation','Country',\
               'PMID','PTyp','Reviews', 'Category','SearchedFor']
  resarray = list()
  paperSummary = Workbook()
  ps = paperSummary.active
  ps.append(['UniprotID', 'Results', 'Synonyms', 'Protein names', 'Category'])
  counter = 2
  for idx in tqdm(range(len(protList))):
    protquer, papercol, papercnt, revcnt = uniProtQuery(protList[idx], IDType, TaxID, \
                                                keywords, KeywordInTitleOnly, \
                                        Mail, maxdate)
    resarray.append([protList[idx], papercnt, protquer['synonyms'], \
                    protquer['protein_names'], protquer['category']])
    ps.append([protList[idx], papercnt, str(protquer['synonyms']), \
                    str(protquer['protein_names']), protquer['category']])
    startcnt = counter+1
    counter+=1
    ps.append(papercol_titles)
    papers = papercol.getPapercolletion()
    for idx in range(papercnt):
      value = list()
      for key in papers.keys():
        value.append(str(papers[key][idx]))
      ps.append(value)
      counter +=1
    ps.row_dimensions.group(startcnt,counter, hidden=True)
    counter += 1

    allRes = pd.DataFrame(resarray, columns = ['UniprotID', 'Results', \
                                                'Synonyms', 'Protein names', \
                                                'Category'])

  if KeywordInTitleOnly:
    paperSummary.save(f'{Filename}_inTitleOnly.xlsx')
  else:
    paperSummary.save(f'{Filename}.xlsx')

  return allRes

def runOPM():
  maxdate= date.today().strftime('%Y/%m/%d')
  papercol_titles = ['Titles','Authors','Journal','Years',\
                             'PMID','DOI','Link']
  resarray = list()
  paperSummary = Workbook()
  ps = paperSummary.active
  ps.append( ['UniprotID', 'Results', 'Reviews', 'Synonyms', 'Protein names', 'Category'])
  counter = 2
  for idx in tqdm(range(len(protList))):
    protquer, papercol, papercnt, revcnt = uniProtQuery(protList[idx], IDType, TaxID, \
                                                keywords, KeywordInTitleOnly, \
                                        Mail, maxdate)
    resarray.append([protList[idx], papercnt, revcnt, protquer['synonyms'], \
                    protquer['protein_names'], protquer['category']])
    ps.append([protList[idx], papercnt, revcnt, str(protquer['synonyms']), \
                    str(protquer['protein_names']), protquer['category']])
    startcnt = counter+1
    counter+=1
    ps.append(papercol_titles)
    papers = papercol.getPapercolletion()
    for idx in range(len(papers['Titles'])):
      value = list()
      for key in papers.keys():
        value.append(str(papers[key][idx]))
      ps.append(value)
      counter +=1
    ps.row_dimensions.group(startcnt,counter, hidden=True)
    counter += 1

    allRes = pd.DataFrame(resarray, columns = ['UniprotID', 'Results', 'Reviews',  \
                                                'Synonyms', 'Protein names', \
                                                'Category'])

  if KeywordInTitleOnly:
    paperSummary.save(f'{Filename}_inTitleOnly.xlsx')
  else:
    paperSummary.save(f'{Filename}.xlsx')

  return allRes

if PubMed == True:
  allRes = runOLM()
else:
  allRes = runOPM()
display(allRes)


pivot_table = allRes.pivot_table(columns=['Category'], aggfunc='size')
color_map = {1: '#ff8000', 2: '#ffc080', 3: '#55a0fb', 4: '#d4d4d4'}
colors = [color_map[cat] for cat in pivot_table.index]
pivot_table.plot.pie(figsize=(6,6),
                     ylabel='',
                     colors=colors)
plt.show()


#Download button
def on_button_clicked(b):
  if KeywordInTitleOnly:
    files.download(f'{Filename}_inTitleOnly.xlsx')
  else:
    files.download(f'{Filename}.xlsx')

button = widgets.Button(description='Download Excel File')
button.on_click(on_button_clicked)
display(button)




