import re
import requests
from requests.adapters import HTTPAdapter, Retry

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


url = 'https://rest.uniprot.org/uniprotkb/search?query=%2A&format=tsv&force=true&fields=accession,protein_name,length,lineage,lineage_ids,organism_name,ft_signal,ft_transmem,xref_tcdb,xref_eggnog,xref_pfam,xref_tigrfams,go_id,xref_interpro,ec,xref_biocyc,ft_dna_bind,ft_binding,cc_subcellular_location,xref_kegg,rhea&compress=true&size=500'
progress = 0
with open('uniprot_download.tsv', 'w') as f:
    for batch, total in get_batch(url):
        lines = batch.text.splitlines()
        if not progress:
            print(lines[0], file=f)
        for line in lines[1:]:
            print(line, file=f)
        progress += len(lines[1:])
        print(f'{progress} / {total}')

