#!/usr/bin/env python
# coding: utf-8

# In[126]:


import requests
from bs4 import BeautifulSoup
import re
  
host = 'https://atmos.uw.edu/~roj/nobackup/MARCUS_and_MICRE/Datasets/MICRE/AAD_radiosondes/'
download_path = '/Volumes/STANFORD_1/micre_soundings/'


req = requests.get(host)
soup = BeautifulSoup(req.text,"lxml")


files = soup.findAll('a',href=re.compile('2016050*'))



type(files[0])



for file in files:
    url = host + file.get('href')
    outfile = url.split('/')
    r = requests.get(url,stream=True)
    with open(download_path+outfile[-1],'wb') as f:
        f.write(r.content)


