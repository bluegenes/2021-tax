#! /usr/bin/env python
"""
Compare "true"/reference annotation with sourmash tax classify
inputs:
  - ref annotation csv: genome,annotation
  - sourmash tax classify csv
"""
import os
import sys
import csv
import argparse
from collections import namedtuple, Counter

import sourmash
from sourmash.lca.lca_utils import (is_lineage_match, taxlist, make_lineage, display_lineage)
from sourmash.tax.tax_utils import (ascending_taxlist, ClassificationResult)
from sourmash.logging import notify, error

CompareClassification = namedtuple('CompareClassification',
                                   'query_name,match_type,match_rank,ref_lineage,ref_type,sourmash_lineage,classification_status,classification_rank,classification_fraction,classification_fname')
CompareReport = namedtuple('CompareReport', 'f_match,f_higher,total_checked,missing_ref,tax_fname')


def get_higher_ranks_to_check(rank):
    remaining_taxlist = []
    tl = list(ascending_taxlist(include_strain=True))
    rank_index = tl.index(rank)
    remaining_taxlist = tl[rank_index+1:]
    return remaining_taxlist


def check_lineages(ref_lin, smash_lin, rank):
    "Check classification for a single genome"
    # match at current the classification rank?
    if is_lineage_match(ref_lin, smash_lin, rank):
        return 'match', rank
    else:
        # match at a higher rank?
        ranks_to_check = get_higher_ranks_to_check(rank)
        for check_rank in ranks_to_check:
            if is_lineage_match(ref_lin,smash_lin,check_rank):
                return 'higher_rank', check_rank
    # no match at all
    return 'nomatch','nomatch'


def load_sourmash_genome_csv(g_csv, *, tax_colnames = ClassificationResult._fields, delimiter=',', force=False):
    '''Load a set of sourmash tax genomes csv_summary results'''
    classif_results = []
    with open(g_csv, 'r') as fp:
        r = csv.DictReader(fp, delimiter=delimiter)
        header = r.fieldnames
        # check for empty file
        if not header:
            raise ValueError(f'Cannot read tax results from {g_csv}. Is file empty?')
        #check for critical column names
        if not set(tax_colnames).issubset(header):
            raise ValueError(f'Did not find all csv_summary column names in {g_csv} - is this an output file from `sourmash tax genome`?')
        for n, row in enumerate(r):
            classif_results.append(row)
    return classif_results


def load_ref_annotations(ref_csv, *, ref_cols = ['genome','lineage'], delimiter=',', refD={}, force=False):
    '''Load a set of reference annotations (genome,lineage)'''
    with open(ref_csv, 'r') as fp:
        r = csv.DictReader(fp, delimiter=delimiter)
        header = r.fieldnames
        # check for empty file
        if not header:
            raise ValueError(f'Cannot read tax results from {g_csv}. Is file empty?')
        #check for correct column names
        if not set(ref_cols).issubset(header):
            raise ValueError(f'Reference annotations ({ref_csv}) must have `genome,lineage` columns to enable proper matching')
        # build reference dictionary
        for n, row in enumerate(r):
            genome = row['genome']
            lin = row['lineage']
            if genome in refD.keys():
                if force:
                    notify(f'Genome {genome} duplicated in provided reference lineages. --force is set, skipping duplicate...')
                    continue
                else:
                    raise ValueError(f'Genome {genome} duplicated provided reference lineages')
            try:
                lin = make_lineage(row['lineage'])
                refD[genome] = lin
            # Not sure if make_lineage will return a ValueError! CHECK.
            except ValueError as exc:
                if force:
                    notify(f'Could not interpret reference lineage row {row}. --force is set, skipping...')
                    continue
                else:
                    raise ValueError(f'Could not interpret reference lineage row {row}. Is the lineage in the right format?')
    return refD


def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    ref_type = args.reference_type
    refD = {}
    # load reference lineages
    for ref_lins in args.ref_annotations:
        try:
            refD = load_ref_annotations(ref_lins, refD=refD, force=args.force)
        except ValueError as exc:
            error(f'ERROR: {str(exc)}')
            sys.exit(-1)

    # keep track of matches
    all_results  = []
    file_reports = []
    all_exact    = 0
    all_higher   = 0
    all_missed   = 0
    all_ref_missed   = 0
    all_cl_missed   = 0
    all_checked  = 0
    higher_counter = Counter()

    # check each tax classify csv
    for n, cl in enumerate(args.tax_genome):
        if n !=0 and n % 10 == 0:
            notify(f"processed {n} classification files.")
        try:
            cl_results = load_sourmash_genome_csv(cl, force=args.force)
        except ValueError as exc:
            error('ERROR: {str(exc)}')
            sys.exit(-1)

        # process each classification result
        num_exact = 0
        num_higher = 0
        ref_missed = 0
        cl_missed = 0
        for m, classif in enumerate(cl_results):
            query = classif['query_name']
            smash_lin = make_lineage(classif['lineage'])
            # get query_lineage
            ref_lin = refD.get(query, None)
            if not ref_lin:
                if args.force:
                    notify(f'No reference lineage provided for query {query}. --force is set, skipping...')
                    match_type='no_ref'
                    match_rank='no_ref'
                    if classif['rank'] == '':
                        match_type='no_ref_or_classif'
                        match_rank='no_ref_or_classif'
                        cl_missed+=1
                    query_result = CompareClassification(query,match_type,match_rank,ref_lin,ref_type,smash_lin,classif["status"],classif['rank'],classif['fraction'],cl)
                    all_results.append(query_result)
                    ref_missed+=1
                    continue
                else:
                    raise ValueError(f'ERROR: No reference lineage provided for query {query}.')
           # handle empty gather results or nomatch
            if classif['rank'] == '':
                notify(f'No classification lineage provided for query {query}. --force is set, skipping...')
                match_type='no_classif'
                match_rank='no_classif'
                query_result = CompareClassification(query,match_type,match_rank,ref_lin,ref_type,smash_lin,classif["status"],classif['rank'],classif['fraction'],cl)
                all_results.append(query_result)
                cl_missed+=1
                continue

            match_type, match_rank = check_lineages(ref_lin, smash_lin, classif['rank'])
            if match_type == 'match':
                num_exact +=1
            elif match_type == 'higher_rank':
                num_higher +=1
                higher_counter[match_rank] +=1
            #namedtuple("ClassificationResult", "query_name, status, rank, fraction, lineage, query_md5, query_filename")
            query_result = CompareClassification(query,match_type,match_rank,ref_lin,ref_type,smash_lin,classif["status"],classif['rank'],classif['fraction'],cl)
            all_results.append(query_result)

        # do some counting / assessment
        num_missed = ref_missed + cl_missed
        num_checked = m - num_missed + 1 # +1 to make one-based???
        fraction_exact = float(num_exact)/num_checked
        fraction_higher = float(num_higher)/num_checked
        fraction_anymatch = float((num_exact+num_higher))/num_checked
        notify(f'for {cl}, {fraction_anymatch} matched reference lineage ({fraction_exact} matched exactly).')

        # store these #
        rep = CompareReport(fraction_exact,fraction_higher,num_checked,num_missed,cl)
        file_reports.append(rep)
        # keep total count too
        all_exact +=num_exact
        all_higher += num_higher
        all_missed +=num_missed
        all_ref_missed +=ref_missed
        all_cl_missed +=cl_missed
        all_checked +=num_checked


    # save query-specific comparison results!
    notify(f"saving lineage comparison {args.output_csv}")
    with open(args.output_csv, 'wt') as fp:
        header = CompareClassification._fields
        w = csv.DictWriter(fp, header, delimiter=',')
        w.writeheader()
        for res in all_results:
            rD = res._asdict()
            rD['classification_fraction'] = f'{float(res.classification_fraction):.3f}'
            if res.ref_lineage:
                rD['ref_lineage'] = display_lineage(res.ref_lineage)
            if res.sourmash_lineage:
                rD['sourmash_lineage'] = display_lineage(res.sourmash_lineage)
            w.writerow(rD)

    ## aggregated reporting
    p_exact = (float(all_exact)/all_checked)*100
    p_higher= (float(all_higher)/all_checked)*100
    p_anymatch = p_exact + p_higher

    notify(f'Compared {all_checked} lineages to {args.reference_type} set.')
    if all_missed:
        notify(f'Note: {all_missed} queries could not be compared.')
        notify(f'Of these, {all_ref_missed} queries did not have a reference lineage; {all_cl_missed} queries did not have a classified lineage.')
    notify(f'\nTotal lineage matches: {p_anymatch:.2f}%')
    notify(f'\tExact Rank: {p_exact:.2f}%')
    notify(f'\tHigher Rank: {p_higher:.2f}%')
    #for rank, count in higher_counter.most_common():
    for rank in ascending_taxlist():
        if rank in higher_counter.keys():
            count = higher_counter[rank]
            p_at_rank = float(count)/all_checked *100
            notify(f'\t  {rank}: {p_at_rank:.2f}%')

    if args.report_csv:
        notify(f"saving per-file comparison report to {args.report_csv}")
        with open(args.report_csv, 'wt') as fp:
            header = CompareReport._fields
            w = csv.DictWriter(fp, header, delimiter=',')
            w.writeheader()
            for rep in file_reports:
                rD = rep._asdict()
                rD['f_match'] = f'{rep.f_match:.3f}'
                rD['f_higher'] = f'{rep.f_higher:.3f}'
                w.writerow(rD)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--tax-genome', required=True, nargs='+') # tax genome results
    p.add_argument('--ref-annotations', required=True, nargs='+')
    p.add_argument('--output-csv', required=True)
    p.add_argument('--report-csv')
    p.add_argument('--reference-type', default = 'reference')
    p.add_argument('--force', action='store_true')
    args = p.parse_args()
    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
