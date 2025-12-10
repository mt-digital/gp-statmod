from datetime import date
from pandas import read_csv, DataFrame

from scholarly import scholarly


ARTICLESDF = read_csv('data/titles.csv')


def make_num_citations_csv(output_file_root = 'data/output/citation_count'):

    outdf = ARTICLESDF.copy()
    outdf['num_citations'] = outdf.Title.map(lambda title: get_num_citations(title))

    outdf.to_csv(output_file_root + str(date.today()) + '.csv', index=False)

    return outdf


def get_num_citations(title):

    with open('data/diagnostic/num_citations_log.txt', 'a') as log_f:
        result = next(scholarly.search_pubs(title))

        # Uncomment to see full result to double check citation count.
        # print(result)

        found_title = result['bib']['title']
        num_citations = result['num_citations']

        # Tell the user what's been found to use as entry for
        diagnostic_msg = \
            f'* Found "{found_title}" result\nfor searched title "{title}"\nwith ' \
            f'{num_citations} citations\n'

        print(diagnostic_msg)

        log_f.write(diagnostic_msg + '\n')

        return num_citations
