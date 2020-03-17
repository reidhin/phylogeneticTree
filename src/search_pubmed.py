# libraries
import os
from Bio import Entrez
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


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


def create_data(data_range: pd.date_range, terms_list: list) -> pd.DataFrame:
    # create data frame
    df_out = pd.DataFrame(
        {
            'mindate': date_range[:-1],
            'maxdate': date_range[1:]
        }
    )

    for term in terms_list:
        print('searching {}'.format(term))
        df_out['{} count'.format(term)] = df_out.apply(
            lambda row: return_count(row['mindate'], row['maxdate'], term),
            axis=1
        )
        df_out['{} cumulative'.format(term)] = pd.to_numeric(df_out['{} count'.format(term)]).cumsum()
    return df_out


if __name__ == '__main__':
    date_range = pd.date_range(start='1920-01-01', end='2020-01-01', freq='AS')
    years = mdates.YearLocator()
    years_fmt = mdates.DateFormatter('%Y')
    data_dir = os.path.join('..', 'data')
    figure_dir = os.path.join('..', 'figures')
    terms = [
        'Multiple Myeloma',
        'Multiple Myeloma AND patient',
        'Multiple Myeloma AND quality of life'
        ]

    if os.path.isfile(os.path.join(data_dir, 'pubmed.csv')):
        print('using local copy')
        df = pd.read_csv(os.path.join(data_dir, 'pubmed.csv'))
        for date_col in ['mindate', 'maxdate']:
            df[date_col] = pd.to_datetime(df[date_col])
    else:
        print('downloading from pubmed')
        df = create_data(date_range, terms)
        df.to_csv(os.path.join(data_dir, 'pubmed.csv'), index=False)

    print(df.head())

    # publications per year
    fig, axs = plt.subplots(1, 1)
    for term in terms:
        df.set_index('mindate')['{} count'.format(term)].plot(
            label=term,
            ax=axs
        )
    axs.legend()
    axs.set_title('Number of yearly publications on PubMed')
    axs.set_ylabel('# / year')
    axs.set_xlabel('Year')
    fig.savefig(os.path.join(figure_dir, 'publications_per_year.png'))

    # cumulative publications
    fig, axs = plt.subplots(2, 1, sharex='all', figsize=(5, 8))
    for term in terms:
        df.set_index('mindate')['{} cumulative'.format(term)].plot(
            label=term,
            ax=axs[0]
        )
    axs[0].legend(title='Search term')
    axs[0].set_title('Number of publications on PubMed')
    axs[0].set_ylabel('#')

    # fraction of publications with QoL
    df['fraction'] = 100 * df['{} count'.format(terms[2])].div(df['{} count'.format(terms[0])], fill_value=0)

    df.set_index('mindate')['fraction'].plot(
        label="Percentage publications on PubMed with\n '{}' ".format(terms[-1]),
        ax=axs[1]
    )
    axs[1].set_ylabel('Percentage with QoL / year')
    axs[1].set_xlabel('Year')
    axs[1].legend()
    plt.tight_layout()
    fig.savefig(os.path.join(figure_dir, 'fraction_per_year.png'), dpi=200)

    plt.show()
