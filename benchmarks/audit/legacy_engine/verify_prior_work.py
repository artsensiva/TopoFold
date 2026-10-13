#!/usr/bin/env python3
"""Verify prior-work citations against Crossref (AGENTS.md rule 9). Records the returned metadata; a reference
is 'verified' only if the returned title matches the expected title (token overlap ≥ 0.8) and the year matches.

Run: ../wpd_local_retrieval/build/venv/bin/python -B verify_prior_work.py
"""
import datetime
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent / 'results' / 'le_prior_work_crossref.json'
UA = {'User-Agent': 'TopoFold-audit/0.1 (citation verification)'}

REFS = [  # key, expected title, year, DOI (None = search by title)
    ('kabsch1976', 'A solution for the best rotation to relate two sets of vectors', 1976, '10.1107/S0567739476001873'),
    ('horn1987', 'Closed-form solution of absolute orientation using unit quaternions', 1987, '10.1364/JOSAA.4.000629'),
    ('theobald2005', 'Rapid calculation of RMSDs using a quaternion-based characteristic polynomial', 2005,
     '10.1107/S0108767305015266'),
    ('levitt1976', 'A simplified representation of protein conformations for rapid simulation of protein folding',
     1976, '10.1016/0022-2836(76)90004-8'),
    ('ramachandran1963', 'Stereochemistry of polypeptide chain configurations', 1963, '10.1016/S0022-2836(63)80023-6'),
    ('oldfield1994', 'Analysis of C alpha geometry in protein structures', 1994, None),
    ('rackovsky1978', 'Differential geometry and polymer conformations. 1. Comparison of protein conformations',
     1978, None),
    ('hu2011', 'Discrete Frenet frame, inflection point solitons, and curve visualization with applications to '
     'folded proteins', 2011, None),
    ('altgodau1995', 'Computing the Fréchet distance between two polygonal curves', 1995, '10.1142/S0218195995000064'),
    ('jiang2008', 'Protein structure-structure alignment with discrete Fréchet distance', 2008, None),
    ('sakoe1978', 'Dynamic programming algorithm optimization for spoken word recognition', 1978,
     '10.1109/TASSP.1978.1163055'),
    ('debrevern2000', 'Bayesian probabilistic approach for predicting backbone structures in terms of protein blocks',
     2000, None),
    ('kolodny2002', 'Small libraries of protein fragments model native protein structures accurately', 2002,
     '10.1016/S0022-2836(02)00942-7'),
    ('vankempen2024', 'Fast and accurate protein structure search with Foldseek', 2024, '10.1038/s41587-023-01773-0'),
    ('mu2005', 'Energy landscape of a small peptide revealed by dihedral angle principal component analysis', 2005,
     '10.1002/prot.20310'),
    ('altis2007', 'Dihedral angle principal component analysis of molecular dynamics simulations', 2007,
     '10.1063/1.2746330'),
    ('rogen2003', 'Automatic classification of protein structure by using Gauss integrals', 2003,
     None),
    ('michaudagrawal2011', 'MDAnalysis: A toolkit for the analysis of molecular dynamics simulations', 2011,
     '10.1002/jcc.21787'),
    ('mcgibbon2015', 'MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories', 2015,
     '10.1016/j.bpj.2015.08.015'),
    ('yianilos1993', 'Data structures and algorithms for nearest neighbor search in general metric spaces', 1993, None),
    ('mirarchi2024', 'mdCATH: A Large-Scale MD Dataset for Data-Driven Computational Biophysics', 2024,
     '10.1038/s41597-024-04140-z'),
]


def tokens(s):
    return set(re.findall(r'[a-z0-9]+', s.lower().replace('é', 'e')))


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
        return json.load(r)['message']


def summarize(m):
    year = None
    for k in ('published-print', 'published-online', 'issued'):
        if m.get(k, {}).get('date-parts'):
            year = m[k]['date-parts'][0][0]
            break
    auth = m.get('author', [{}])
    return {'doi': m.get('DOI'), 'title': (m.get('title') or [''])[0], 'year': year,
            'first_author': (auth[0].get('family') if auth else None),
            'container': (m.get('container-title') or [''])[0], 'volume': m.get('volume'), 'page': m.get('page')}


def main():
    out = {'_meta': {'checked_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'source': 'api.crossref.org',
                     'rule': 'verified if title token overlap >= 0.8 and year within +-1 (print vs online)'},
           'refs': {}}
    for key, title, year, doi in REFS:
        rec = {'expected_title': title, 'expected_year': year, 'queried_doi': doi}
        try:
            if doi:
                m = summarize(get('https://api.crossref.org/works/' + urllib.parse.quote(doi)))
            else:
                q = urllib.parse.urlencode({'query.bibliographic': title, 'rows': 1})
                res = json.load(urllib.request.urlopen(urllib.request.Request(
                    'https://api.crossref.org/works?' + q, headers=UA), timeout=60))['message']['items']
                m = summarize(res[0]) if res else None
            rec['crossref'] = m
            if m:
                a, b = tokens(title), tokens(m['title'])
                overlap = len(a & b) / max(len(a), 1)
                rec['title_overlap'] = round(overlap, 3)
                rec['verified'] = bool(overlap >= 0.8 and m['year'] is not None and abs(m['year'] - year) <= 1)
            else:
                rec['verified'] = False
        except Exception as e:  # noqa: BLE001
            rec['error'] = repr(e)
            rec['verified'] = False
        out['refs'][key] = rec
        print(key, rec['verified'], rec.get('crossref', {}) and rec['crossref'].get('doi'), rec.get('title_overlap'))
        time.sleep(1.5)
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')


if __name__ == '__main__':
    main()
