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
from openpyxl import Workbook
import ipywidgets as widgets
from google.colab import files
import nltk
nltk.download('punkt')

import copy
#@title Install dependencies and enter input
#@markdown  Please enter your E-Mail (Needed for the search, you will not receive any E-Mails)
Mail = 'an.gocke@uke.de'#@param {type:"string"}
#@markdown  Specify the title of the output file. If you have searched in the title only, it will be added to the filename
Filename = 'Suppl10_OLM2_2024'#@param {type:"string"}

Protein_List = 'THBS2 CAV2 SCG2 SLC6A1 SAV1 SEZ6L2 ERO1A RAB3B OBSL1 CD109 PTPN14 MRPL35 LRPAP1' #@param {type:"string"}
#@markdown  Specify the type of the Input for the Protein List, set a tick either for the Accession or the Gene Name
Accession = False #@param {type:"boolean"}
Gene_name = True #@param {type:"boolean"}


if (Accession and Gene_name) or (not Accession and not Gene_name):
  raise RuntimeError(f'You need to enter either Accession or Gene_name')

#@markdown  Taxonomy ID (e.g. homo sapiens: 9606; mus musculus: 10090)
TaxID = "9606"#@param {type:"string"}

#@markdown  If you want to enter multiple keywords, please separate them by semicolons. If multiple keywords should be found together in the abstract or title, please tick the "Together" button. Else the tool will search whether one of the keywords was found.
Together = False #@param {type:"boolean"}
Keywords = 'migration'#@param {type:"string"}

#@markdown  Checks if the keywords are mentioned  in the title only or in the title and the abstract
KeywordInTitleOnly = False #@param {type:"boolean"}

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
  start = "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:"
  middle = '+AND+accession:'\
          if str(idtype).lower() == 'accession' \
          else "+AND+gene:"
  end = '&format=tsv'
  contents = urllib.request.urlopen('{}{}{}{}{}'.format(start,
                taxid, middle, uniprotID, end)).read().decode("utf-8")
  data = io.StringIO(contents)
  df = pd.read_csv(data, sep="\t", header=0)
  return df

def getUniProtSynonyms(uniprotID, idtype, taxid):
  contents = requestOnUniProtWebpage(uniprotID, idtype, taxid, False)
  if len(contents['Gene Names']) > 0:
    syns = [str(x) for x in contents['Gene Names']][0]
    syns = re.sub(r'\s', ',', syns)
    syns = syns.split(",")
  else:
    syns = ''
  if len(contents['Protein names']) > 0:
    prots = [str(x) for x in contents['Protein names']][0]
    contents = re.split(r'\)?\s[\(\[]', re.sub(r'\]|\)', '', prots))
    for idx in range(len(contents)):
      cur_str = re.sub(r'^\s', '', contents[idx])
      cur_str = re.sub(r'\s', '-', cur_str)
      contents[idx] = re.sub(r'\"|\;', '', cur_str)
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
  searchsummary['keywordInTitleOnly'] = tiab
  searchsummary['totalResults'] = 0
  searchsummary['category'] = 0
  searchsummary['false'] = 0
  papercol = Papercollection(uniprotID, idtype, taxid, synonyms, protnames,\
                             mail, KeywordInTitleOnly, keyword,  maxdate)
  searchsummary['category'] = papercol.bestCategory()
  searchsummary['totalResults'] = papercol.resultcnt()
  return searchsummary, papercol, papercol.resultcnt()


class Papercollection:
  def __init__(self, uniprotID, idtype, taxid, synonyms, protnames,
                mail, titleonly, keyword, maxdate):
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
    self._bestcategory = 4
    self._PMIDdict = dict()
    self._maxdate = maxdate
    self._keywordlist = keyword.split(';')

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

    if self._bestcategory == 4:
      if len(self._syns) < 1 and len(self._prots) < 1:
        self._bestcategory = 0
      else:
        self._bestcategory = 3
      if  self.resultcnt() != 0:
        self._bestcategory = min(self._papercollection['Category'])


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

  def bestCategory(self):
    return self._bestcategory

  def resultcnt(self):
    if self._papercollection['Titles'] == ['nan']:
      x = 0
    else:
      x = len(self._papercollection['Titles'])
    return x

  def fetch_handle(self, id_list):
    Entrez.email = self._mail
    handle = Entrez.efetch(db='pubmed', id=id_list, rettype='medline',
                                                    retmode='text')
    return(Medline.parse(handle))

  def searchAll(self,synonyms, proteinnames):
    syns = []
    if len(synonyms) > 0 or len(proteinnames) > 0:
      syns =  synonyms + proteinnames
    syns.append(self._uniprotID)
    time.sleep(0.3)
    if len(syns) > 0:
      results, term = self.searchAllSyns(syns)
      id_list = results['IdList']
    return(self.fetch_handle(id_list), term)

  def searchAllSyns(self, syns):
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


def runMain():
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
    protquer, papercol, papercnt = uniProtQuery(protList[idx], IDType, TaxID, \
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

allRes = runMain()
display(allRes)


pivot_table = allRes.pivot_table(columns=['Category'], aggfunc='size')
color_map = {0: '#d4d4d4', 1: '#ff8000', 2: '#ffc080', 3: '#55a0fb'}
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




