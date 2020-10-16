import requests
from bs4 import BeautifulSoup
import re
  
host = 'https://atmos.uw.edu/~roj/nobackup/MARCUS_and_MICRE/Datasets/MICRE/Cloud_and_Precipitation_Properties/V1.6/figures/'
download_path = '/Volumes/STANFORD_6/cloud_precip_properties/'

req = requests.get(host)
soup = BeautifulSoup(req.text,"lxml")

files = soup.findAll('a',href=re.compile('MICRE_AAD_'))
#files = soup.findAll('a')
#print(files[0])

for file in files:
    print(file)
    url = host + file.get('href')
    outfile = url.split('/')
    r = requests.get(url,stream=True)
    with open(download_path+outfile[-1],'wb') as f:
        f.write(r.content)






