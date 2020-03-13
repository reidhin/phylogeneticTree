# libraries
from Bio import Entrez
import pandas as pd
import matplotlib.pyplot as plt


def return_count(mindate, maxdate, term):
    Entrez.email = 'hans@orikami.nl'
    h = Entrez.esearch(
        db='pubmed',
        term=term,
        mindate=mindate.strftime('%Y/%m/%d'),
        maxdate=maxdate.strftime('%Y/%m/%d'),
        datetype='pdat',
        rettype='count'
    )
    result = Entrez.read(h)
    print(result['Count'])
    return result['Count']


date_range = pd.date_range(start='1915-01-01', end='2020-01-01', freq='AS')
# create data frame
df = pd.DataFrame(
    {
        'mindate': date_range[:-1],
        'maxdate': date_range[1:]
    }
)

TERM = 'Multiple Myeloma'

df['Count'] = df.apply(lambda row: return_count(row['mindate'], row['maxdate'], TERM), axis=1)

df['cumsum'] = pd.to_numeric(df['Count']).cumsum()

print(df.head())

df.plot(x='mindate', y='cumsum')
plt.show()
