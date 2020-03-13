# libraries
from Bio import Entrez
import pandas as pd
import matplotlib.pyplot as plt


def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='20',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results


def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'hans@orikami.nl'
    handle = Entrez.efetch(db='pubmed',
                           retmode='medline',
                           id=ids)
    results = Entrez.read(handle)
    return results


if __name__ == '__main__':
    from Bio import Entrez
    from Bio import Medline

    MAX_COUNT = 2000
    TERM = 'Multiple Myeloma'

    print('Getting {0} publications containing {1}...'.format(MAX_COUNT, TERM))
    Entrez.email = 'hans@orikami.nl'
    h = Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=TERM)
    result = Entrez.read(h)
    print('Total number of publications containing {0}: {1}'.format(TERM, result['Count']))
    ids = result['IdList']
    h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
    records = Medline.parse(h)

    titles, DP, DEP, EDAT = [], [], [], []
    for record in records:
        titles.append(record.get('TI', ''))
        DP.append(record.get('DP', ''))
        DEP.append(record.get('DEP', ''))
        EDAT.append(record.get('EDAT', ''))

    df = pd.DataFrame({
        'Title': titles,
        'DP': pd.to_datetime(DP, format='%Y %b %d', errors='coerce'),
        'DEP': pd.to_datetime(DEP, format='%Y%m%d', errors='coerce'),
        'EDAT': pd.to_datetime(EDAT, format='%Y/%m/%d %H:%M', errors='coerce')
    })
    df = df.sort_values(by='EDAT')
    df['cumsum'] = range(1, len(df) + 1)
    df.plot(x='EDAT', y='cumsum')
    df.plot(x='DEP', y='cumsum')

    print(df.head())
    '''
    results = search('Multiple Myeloma')
    id_list = results['IdList']
    papers = fetch_details(id_list)
    titles, years, months, days = list(), list(), list(), list()
    for i, paper in enumerate(papers['PubmedArticle']):
        print("%d) %s" % (i+1, paper['MedlineCitation']['Article']['ArticleTitle']))
        titles = titles.append(paper['MedlineCitation']['Article']['ArticleTitle'])
    # Pretty print the first paper in full to observe its structure
    #import json
    #print(json.dumps(papers[0], indent=2, separators=(',', ':')))
    '''