"""
Parse ClinVar vcf and append extra information from the ClinVar XML
"""
import argparse
import gzip
import io
import logging
import os
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from csv import writer
from datetime import datetime
import pandas as pd


def expand_clinvar_vcf(xml_file, vcf_file, out_file, pre_may_2017=False):
    """
    Args:
        xml_file: input ClinVar XML annotation file
        vcf_file: input ClinVar VCF file
        out_file: output annotated VCF
        pre_may_2017: assume old format ClinVar file
    """

    logging.info('Reading in vcf file')
    header, vcf_df = read_vcf(vcf_file)

    logging.info('Writing out updated vcf header')

    write_vcf_header(header, out_file)

    logging.info('Creating dictionary for looking up IDs')
    if pre_may_2017:
        logging.info('VCF dataframe size:%s', vcf_df.shape)
        logging.info('Expanding multiallelic variants')
        vcf_df = split_multi_vcf(vcf_df)
        logging.info('VCF dataframe size post expansion:%s', vcf_df.shape)
        # (Only works when RCV IDs are in INFO vcf column - field CLNACC.
        # Works with ClinVar vcf older than May 2017)
        id_dict = split_vcf_info(vcf_df)
    else:
        # (Only works when ClinVar Variation IDs are in ID vcf column.
        # ClinVar vcf older than May 2017 won't work)
        id_dict = dict(zip(vcf_df['ID'], vcf_df.index))

    logging.info('Mining through XML file')
    xml_dict = defaultdict(lambda: defaultdict(list))
    n = 0

    context = ET.iterparse(get_handle(xml_file), events=('start', 'end'))
    for event, elem in context:
        if elem.tag != 'ClinVarSet' or event != 'end':
            continue
        ms_id = elem.find('.//MeasureSet').attrib.get('ID')
        key_id = (elem.find('.//ClinVarAccession').attrib.get('Acc')
                  if pre_may_2017 else ms_id)
        cva = elem.findall('./ClinVarAssertion')

        if key_id in id_dict:
            xml_dict[id_dict[key_id]]['CLNSUBA'] += get_submitters(elem)
            xml_dict[id_dict[key_id]]['CLNREVSTATA'] += get_clinsig(
                    cva, field='status_ordered', record_id=ms_id)
            xml_dict[id_dict[key_id]]['CLNDATEA'] += get_clinsig(
                    cva, field='last_eval', record_id=ms_id)
            xml_dict[id_dict[key_id]]['CLNDATESUBA'] += get_submitdate(
                    cva, field='submitterDate', record_id=ms_id)
            xml_dict[id_dict[key_id]]['CLNSIGA'] += get_clinsig(
                    cva, field='description', record_id=ms_id)
            xml_dict[id_dict[key_id]]['CLNORA'] += get_origin(
                    cva, as_set=False)
            xml_dict[id_dict[key_id]]['CLNCOMA'] += get_clinsig(
                    cva, field='comment', record_id=ms_id)

            # Duplicate disease name if multiple SCV entries in one ClinVarSet
            # record
            scv_ids = get_accession(cva, field='SCV')
            xml_dict[id_dict[key_id]]['CLNSCVA'] += scv_ids
            disease_name = '/'.join(get_traits(elem))
            for _ in range(len(scv_ids)):
                xml_dict[id_dict[key_id]]['CLNDNA'] += [disease_name]

        elem.clear()

        # Progress report
        n += 1
        if n % 50000 == 0:
            logging.info('%s records processed', n)

    logging.info('Finished mining through XML file')
    logging.info('Processed %s ClinVarSet records', n)

    # Add extra column to vcf dataframe with info extracted from xml
    for key, value in xml_dict.items():
        xml_dict[key] = join_entries(value)  # join xml entries into strings
        vcf_df.at[key, 'INFO_add'] = ';'.join(
            [k + '=' + v for k, v in xml_dict[key].items()])

    logging.info('Combining original INFO column with info extracted from xml')
    vcf_df['INFO_updated'] = vcf_df['INFO'] + ';' + vcf_df['INFO_add']
    vcf_df['INFO_updated'] = vcf_df['INFO_updated'].fillna(vcf_df['INFO'])
    vcf_df = vcf_df.drop(columns=['INFO_add', 'INFO']).rename(
        columns={'INFO_updated': 'INFO'})

    logging.info('Writing out the main VCF body to output file')
    # using append mode because header already written
    vcf_df.to_csv(out_file, mode='a', index=False, sep='\t', header=True)


def get_accession(elem, field='SCV', record_id=None):
    """
    Extracts SCV['Acc'] from ClinVarAccession elements
    """
    acc_list = [x.find('./ClinVarAccession') for x in elem]
    val_list = []
    try: # FIXME
        for acc in acc_list:
            if acc.attrib.get('Type') == field:
                val_list.append(acc.attrib.get('Acc'))
    except AttributeError:
        logging.warning('Record %s has no %s record', record_id, field)
    return val_list


def get_clinsig(elem, field=None, count=False, record_id=None):
    """
    Args:
        elem: parent XML element
        field: specifies the element to interrogate - 'status_ordered', 'description','last_eval','comment'
        count: optionally count description
        record_id: supply for detailed feeback in logging

    """
    clinsig_list = [item for x in elem for item in x.findall('./ClinicalSignificance')]
    
    if field == 'status_ordered':
        status_ordered = []
        if len(clinsig_list) == 0:
            logging.warning('Record %s has ClinicalSignificance in %s',
                            record_id, elem[0].tag)
        else:
            try:
                status_ordered = [x.find('./ReviewStatus').text for x in clinsig_list]
            except AttributeError:
                logging.warning('Record %s has no ReviewStatus field in %s',
                                record_id, elem[0].tag)
        return status_ordered

    if field == 'description':
        desc_ordered = []
        if len(clinsig_list) == 0:
            logging.warning('Record %s has no ClinicalSignificance in %s',
                            record_id, elem[0].tag)
        else:
            try:
                desc_ordered = [x.find('./Description').text for x in clinsig_list]
            except AttributeError:
                logging.warning('Record %s has no Description field in %s',
                                record_id, elem[0].tag)
        if count:
            desc_ordered_low = [x.lower() for x in desc_ordered]  # needs to be lower case
            desc_count = {}
            for k in ['pathogenic',
                      'likely_pathogenic',
                      'uncertain_significance',
                      'benign',
                      'likely_benign']:
                k_count = k.replace('_', ' ')  # key for count
                desc_count[k] = str(desc_ordered_low.count(k_count))
            return desc_count
        return desc_ordered

    if field == 'last_eval':
        date_ordered = []
        if len(clinsig_list) == 0:
            logging.warning('Record %s has no ClinicalSignificance in %s',
                            record_id, elem[0].tag)
        else:
            for clinsig in clinsig_list:
                date_ordered.append(clinsig.attrib.get('DateLastEvaluated', '0000-00-00'))
        return date_ordered

    if field == 'comment':
        comment_ordered = []
        if len(clinsig_list) == 0:
            logging.warning('Record %s has no ClinicalSignificance in %s',
                            record_id, elem[0].tag)
        else:
            for clinsig in clinsig_list:
                try:
                    comment_ordered.append(clinsig.find('./Comment').text)
                except AttributeError:
                    comment_ordered.append('')
            return comment_ordered

    raise NotImplementedError(f'field {field} is not handled')


def get_handle(fname, ftype='xml'):
    """
    Return an opened file handle

    Args:
        fname: name of the file
        ftype: 'xml' or 'vcf'
    """
    if fname[-3:] == '.gz' and ftype == 'xml':
        handle = gzip.open(fname)
    elif fname[-3:] == '.gz' and ftype == 'vcf':
        handle = gzip.open(fname, 'rt')
    else:
        handle = open(fname)
    return handle


def get_origin(elem, as_set=True):
    """
    Extract ObservedIn/Sample/Origin
    """
    origin_list = []
    for node in elem:
        try:
            origin_list.append(node.find('./ObservedIn/Sample/Origin').text)
        except AttributeError:
            origin_list.append('')
    if as_set:
        return set(origin_list)
    return origin_list


def get_submitdate(elem, field='submitterDate', record_id=None):
    """
    Extracts submitterDate from ClinVarSubmissionID
    """
    subid_list = [x.find('./ClinVarSubmissionID') for x in elem]
    val_list = []
    if len(subid_list) == 0:
        logging.warning('Record %s has no ClinVarSubmissionID field',
                        record_id)
    else:
        for subid in subid_list:
            val_list.append(subid.attrib.get(field, '0000-00-00'))
    return val_list


def get_submitters(elem, ordered=True):
    """
    Return the set of submitter attribute values from ClinVarSubmissionID
    elements
    """
    submitters = []
    for submitter_node in elem.findall('.//ClinVarSubmissionID'):
        try:
            submitters.append(submitter_node.attrib['submitter'])
        except KeyError:
            return None
    if ordered:
        return submitters
    return set(submitters)


def get_traits(elem):
    """
    Interrogate TraitSet elements for 'traits'
    """
    # now find the disease(s) this variant is associated with
    trait_values = []
    for traitset in elem.findall('.//TraitSet'): 
        disease_name_nodes = traitset.findall('.//Name/ElementValue')
        trait_values.extend([x.text for x in disease_name_nodes
                             if x.attrib.get('Type') == 'Preferred'])
        return trait_values


def join_entries(dct, with_sep='|'):
    """
    Args:
        dct: dictionary
    Returns:
        dct with values updated in-place with joined entries
    """
    for key, value in dct.items():
        if isinstance(value, set):
            value = sorted(value)
        if isinstance(value, (set, list)):
            value = [
                replace_sep(string=v, sep=[';', '|', ' ', '"', ','],
                            replace_with=[':', ':', '_', '', '&'])
                for v in value
            ]
            dct[key] = remove_newlines_and_tabs(with_sep.join(value))
        if value is None:
            dct[key] = ''
    return dct


def read_vcf(fname):
    """
    Read a VCF file.

    Args:
        fname: name of the file
    Returns:
        (header, vcf) where header is a list of header lines and vcf is a
            pandas dataframe
    """
    header = []
    lines = []
    with get_handle(fname, ftype='vcf') as f:
        for line in f:
            if line.startswith('##'):
                header.append(line.strip())
            else:
                lines.append(line)
    vcf = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    )
    return (header, vcf)


def remove_newlines_and_tabs(string):
    """
    Replace each tab and newline character with a single space
    """
    return re.sub('[\t\n\r]', ' ', string)


def replace_sep(string, sep, replace_with):
    """
    Args:
        string: the target string to perform replacement on
        sep: list of separators to replace
        replace_with: list of separators to use as replacement
    Returns:
        string with replacements made
    """
    
    for i, sep_0 in enumerate(sep):
        string = string.replace(sep_0, replace_with[i])

    return string


def split_multi_vcf(df):
    """
    Expand multiallelic variants in old-format VCF dataframe (--pre-may-2017)

    Args:
        df: old-format VCF dataframe

    Returns:
        dataframe with one record per ALT
    """
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)  # change to None for new pandas
    # Write temp csv for splitting multiallelic variants into multiple rows
    # (to increase speed)
    temp_output = io.StringIO()
    csv_writer = writer(temp_output)

    for i, alts in enumerate(df.get('ALT')):
        temp_row = df.iloc[i].copy()
        if ',' in alts:
            for n, alt in enumerate(alts.split(','), 1):
                temp_row['ALT'] = alt
                new_info = []
                for field in df.iloc[i]['INFO'].split(';'):  # extract INFO field for each line
                    if field.startswith('CLNALLE='):
                        alt_info_n = [
                            int(x)
                            for x in field.replace('CLNALLE=', '').split(',')]
                        new_info.append('CLNALLE=1')
                    elif field.startswith('CLN'):
                        field_id, values = field.split('=')
                        if n in alt_info_n:
                            alt_info_index = alt_info_n.index(n)
                            value = values.split(',')[alt_info_index]
                            new_info.append(f'{field_id}={value}')
                        else:
                            new_info.append(f'{field_id}=.')
                    else:
                        new_info.append(field)
                temp_row['INFO'] = ';'.join(new_info)
                csv_writer.writerow(temp_row)
        else:
            csv_writer.writerow(temp_row)

    temp_output.seek(0)  # we need to get back to the start of the StringIO
    df_split = pd.read_csv(temp_output, header=None, names=df.columns)
    return df_split


def split_vcf_info(df):
    """
    Returns a dictionary of {id: line_number}

    Args:
        df: dataframe with INFO column
    """
    info_dct = {}
    for i, info in enumerate(df.get('INFO')):
        for field in info.split(';'):
            if field.startswith('CLNACC='):
                for rcv_id_v in field.split('=')[1].split('|'):  # can have multiple IDs
                    rcv_id = rcv_id_v.split('.')[0]  # removing version
                    if not rcv_id.startswith('RCV') and rcv_id != '':
                        logging.warning(
                        'VCF %s at index %s is not correctly formatted, '
                        'should be CLNACC=RCVXXX', field, i)
                    info_dct[rcv_id] = i  # dictionary with RCV ids and line numbers
    return info_dct


def write_vcf_header(header, out_file):
    """
    Write VCF header
    """
    vcf_header_dict = {
        'CLNSUBA': 'Submitters - all, ordered',
        'CLNSIGA': 'Clinical significance - all, ordered',
        'CLNDATEA': 'Date pathogenicity last reviewed - all, ordered',
        'CLNDATESUBA': 'Submission date - all, ordered',
        'CLNREVSTATA': 'ClinVar review status for the Variation ID - all, ordered',
        'CLNORA': 'Allele origin - all, ordered',
        'CLNSCVA': 'SCV IDs - all, ordered',
        'CLNDNA': 'Preferred disease name - all, ordered',
        'CLNCOMA': 'Comment on clinical significance - all, ordered'
    }

    for key, value in vcf_header_dict.items():
        add_line = f'##INFO=<ID={key},Number=.,Type=String,Description="{value}">'
        header.append(add_line)

    with open(out_file, 'w') as fout:
        fout.write('\n'.join(header))
        fout.write('\n')


def main():
    """
    Parse the command-line and do the work.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    req_grp = parser.add_argument_group(title='Required')
    req_grp.add_argument('-x',
                         '--xml',
                         type=str,
                         help='ClinVar XML file (can be .gz)',
                         required=True)
    req_grp.add_argument('-i',
                         '--input',
                         type=str,
                         help='ClinVar input vcf file name (can be .gz)',
                         required=True)
    req_grp.add_argument('-o',
                         '--out',
                         type=str,
                         help='Output vcf file name (non gz)',
                         required=True)
    parser.add_argument('-l', '--log', type=str, help='Log file name')
    parser.add_argument('--pre-may-2017',
                        action='store_true',
                        default=False,
                        help='Assume old ClinVar file format: '
                             'RCV IDs are in INFO vcf column - field CLNACC')

    args = parser.parse_args()

    if args.log:
        logging.basicConfig(filename=args.log, level=logging.INFO)
    else:
        out_path = os.path.dirname(args.out)
        logging.basicConfig(filename=(out_path + '/clinvar_vcf_parser.log'),
                            level=logging.WARNING)

    logging.info('Start time: %s\n', datetime.now())
    logging.info(args)

    expand_clinvar_vcf(xml_file=args.xml, out_file=args.out,
                       vcf_file=args.input, pre_may_2017=args.pre_may_2017)

    logging.info('\nEnd time: %s', datetime.now())


if __name__ == '__main__':
    main()
