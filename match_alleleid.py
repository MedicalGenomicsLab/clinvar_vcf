import re
import sys
import io
import os
import gzip
import argparse
from collections import defaultdict
import xml.etree.ElementTree as ET
import functools
import operator
import pandas as pd
from datetime import datetime
import logging


def get_handle(file, type ='xml'):
    if file[-3:] == '.gz' and type =='xml':
        handle = gzip.open(file)
    elif file[-3:] == '.gz' and type =='vcf':
        handle = gzip.open(file, 'rt')
    else:
        handle = open(file)
    return handle

def read_vcf(file):
    header =[]
    lines = []
    with get_handle(file, type='vcf') as f:
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
    return (header,vcf)

def reduce_list(l):
    return functools.reduce(operator.concat, l) 

def replace_sep(string, sep=';', replace_with=":"):
    # sep and replace_with can both be lists or str
    if type(sep) == list and type(replace_with) == list:
        for i, sep_0 in enumerate(sep):
            string = string.replace(sep_0, replace_with[i])
    elif type(sep) == list and type(replace_with) == str:
        for sep_0 in sep:
            string = string.replace(sep_0, replace_with)  
    elif type(sep) == str and type(replace_with) == list:
        string = string.replace(sep, replace_with[0])
    else:
        string = string.replace(sep, replace_with)
    return string

def remove_newlines_and_tabs(string):
    return re.sub("[\t\n\r]", " ", string)

def join_entries(d, with_sep='|'):
    for key,value in d.items():
        if type(value) == set:
            value = sorted(value)
        if type(value) == set or type(value) == list:
            value = [replace_sep(v, sep=[';','|',' ','"'], replace_with=[':',':','_','']) for v in value]
            d[key]=remove_newlines_and_tabs(with_sep.join(value))
        if value is None:
            d[key]=''
    return d

def get_submitters(elem, ordered=True):
    submitters= []
    for submitter_node in elem.findall('.//ClinVarSubmissionID'):
        try:
            submitters.append(submitter_node.attrib['submitter'])
        except ValueError:
            return None
    if ordered:
        return submitters 
    else:
        return set(submitters)

def get_clinsig(elem, field=None, count=False, record_id=None):
    clinsig_list = reduce_list([x.findall('./ClinicalSignificance') for x in elem])
    if field == "status": 
        review_status = ''
        if len(clinsig_list) == 0:
            logging.warning(f"Record {record_id} has no ClinicalSignificance in {elem[0].tag}")
        else:
            try:
                review_status = clinsig_list[0].find('./ReviewStatus').text
            except:
                logging.warning(f"Record {record_id} has no ReviewStatus field in {elem[0].tag}")                
        return review_status
    elif field == "status_ordered":
        status_ordered = []
        if len(clinsig_list) == 0:
            logging.warning(f"Record {record_id} has ClinicalSignificance in {elem[0].tag}")
        else:
            try:
                status_ordered = [x.find('./ReviewStatus').text for x in clinsig_list]
            except:
                logging.warning(f"Record {record_id} has no ReviewStatus field in {elem[0].tag}")                  
        return status_ordered
    elif field == "description":
        desc_ordered = []
        if len(clinsig_list) == 0:
            logging.warning(f"Record {record_id} has no ClinicalSignificance in {elem[0].tag}")
        else:
            try:   
                desc_ordered = [x.find('./Description').text for x in clinsig_list]
            except:
                logging.warning(f"Record {record_id} has no Description field in {elem[0].tag}")  
        if count:
            desc_ordered_low = [x.lower() for x in desc_ordered] #needs to be lower case
            desc_count = {}
            for k in ['pathogenic','likely_pathogenic','uncertain_significance','benign','likely_benign']:
                k_count = k.replace('_',' ') #key for count
                desc_count[k] = str(desc_ordered_low.count(k_count))
            return desc_count
        elif not count:
            return desc_ordered
    elif field == "last_eval":
        if len(clinsig_list) == 0:
            logging.warning(f"Record {record_id} has no ClinicalSignificance in {elem[0].tag}")
        else:            
            date_ordered = [x.attrib.get('DateLastEvaluated', '0000-00-00') for x in clinsig_list]
            if len(date_ordered) == 0:
                return '0000-00-00'
            else:
                return date_ordered
    elif field == "comment":
        comment_ordered = []
        if len(clinsig_list) == 0:
            logging.warning(f"Record {record_id} has no ClinicalSignificance in {elem[0].tag}")
        else:
            for clinsig in clinsig_list:
                try:
                    comment_ordered.append(clinsig.find('./Comment').text)
                except:
                    comment_ordered.append('')
            return comment_ordered
         
def get_accession(elem, field=None, key='Acc', record_id=None):
    acc_list = [x.find('./ClinVarAccession') for x in elem]
    try: 
        val_list = []
        for acc in acc_list:
            if acc.attrib.get('Type') == field:
                val_list.append(acc.attrib.get(key))
        return val_list
    except:
        logging.warning(f"Record {record_id} has no {field} record")

def get_traits(elem, field=None):
    # now find the disease(s) this variant is associated with
    for traitset in elem.findall('.//TraitSet'):
        if field == "traits":
            trait_values = []
            disease_name_nodes = traitset.findall('.//Name/ElementValue')
            trait_values = [x.text for x in disease_name_nodes if x.attrib.get('Type') == 'Preferred']
            return trait_values
        elif field == "mechanism":
            attribute_nodes = traitset.findall('.//AttributeSet/Attribute')
            disease_mechanism = {x.text.strip() for x in attribute_nodes if x.attrib.get('Type') == 'disease mechanism'}
            return(disease_mechanism)
        elif field == "xrefs":
            # put all the cross references one column, it may contains NCBI gene ID, conditions ID in disease databases.
            xref_nodes = traitset.findall('.//XRef')
            xrefs = {f"{x.attrib.get('DB')}:{x.attrib.get('ID')}".replace(';',',') for x in xref_nodes}
            return(xrefs)
                     
def get_origin(elem, as_set=True):
    origin_list=[]
    for node in elem:
        try:
            origin_list.append(node.find('./ObservedIn/Sample/Origin').text)
        except:
            origin_list.append('')
    if as_set:                 
        return set(origin_list)
    else:
        return origin_list

def write_vcf_header(header, out_file):
    vcf_header_dict={'CLNSUBA':'Submitters - all, ordered',
                'CLNSIGA':'Clinical significance - all, ordered',
                'CLNDATEA':'Date pathogenicity last reviewed - all, ordered',
                'CLNDATEUPA':'Date accession last updated - all, ordered',
                'CLNREVSTATA':'ClinVar review status for the Variation ID - all, ordered',
                'CLNORA': 'Allele origin - all, ordered',
                'CLNSCVA':'SCV IDs - all, ordered',
                'CLNDNA': 'Preferred disease name - all, ordered',
                'CLNCOMA': 'Comment on clinical significance - all, ordered'}
    
    for key, value in vcf_header_dict.items():
        add_line = f"##INFO=<ID={key},Number=.,Type=String,Description=\"{value}\">"
        header.append(add_line)
    
    with open(out_file,'w') as fout:
        fout.write('\n'.join(header))
        fout.write('\n')
      
def expand_clinvar_vcf(xml_file, vcf_file, out_file):
    
    logging.info('Reading in vcf file')
    header,vcf_df = read_vcf(vcf_file)
    
    logging.info('Writing out updated vcf header')
    
    write_vcf_header(header,out_file)
    
    logging.info('Creating dictionary for looking up IDs')
    # (Only works when ClinVar Variation IDs are in ID vcf column. ClinVar vcf older than May 2017 won't work)
    id_dict=dict(zip(vcf_df['ID'],vcf_df.index))
     
    logging.info('Mining through XML file')
    xml_dict=defaultdict(lambda: defaultdict(list))
    n=0
    
    context = ET.iterparse(get_handle(xml_file), events=("start", "end"))
    for event, elem in context:
        if elem.tag != 'ClinVarSet' or event != 'end':
            continue 
        ms_id = elem.find('.//MeasureSet').attrib.get('ID') 
        cva=elem.findall("./ClinVarAssertion")

        if ms_id in id_dict:
            xml_dict[id_dict[ms_id]]['CLNSUBA'] += get_submitters(elem)
            xml_dict[id_dict[ms_id]]['CLNSIGA'] += get_clinsig(cva, field='status_ordered', record_id=ms_id)
            xml_dict[id_dict[ms_id]]['CLNDATEA'] += get_clinsig(cva, field='last_eval', record_id=ms_id) 
            xml_dict[id_dict[ms_id]]['CLNDATEUPA'] += get_accession(cva, field='SCV', key='DateUpdated')
            xml_dict[id_dict[ms_id]]['CLNREVSTATA'] += get_clinsig(cva, field='description', record_id=ms_id)        
            xml_dict[id_dict[ms_id]]['CLNORA'] += get_origin(cva, as_set=False) 
            xml_dict[id_dict[ms_id]]['CLNCOMA'] += get_clinsig(cva, field='comment', record_id=ms_id)  
            
            # Duplicate disease name if multiple SCV entries in one ClinVarSet record
            scv_ids = get_accession(cva, field='SCV')
            xml_dict[id_dict[ms_id]]['CLNSCVA'] += scv_ids  
            disease_name = '/'.join(get_traits(elem, field='traits'))
            for _ in range(len(scv_ids)):
                xml_dict[id_dict[ms_id]]['CLNDNA'] += [disease_name]
 
        elem.clear()
        
        #Progress report
        n+=1
        if n % 50000 == 0:
            logging.info(f"{n} records processed") 

    logging.info('Finished mining through XML file')
    logging.info(f"Processed {n} ClinVarSet records")
    
    
    # Add extra column to vcf dataframe with info extracted from xml
    for key,value in xml_dict.items():
        xml_dict[key] = join_entries(value) # join xml entries into strings 
        vcf_df.at[key, 'INFO_add'] = ';'.join([k + "=" + v for k, v in xml_dict[key].items()])
    
    logging.info('Combining original INFO column with info extracted from xml')
    vcf_df['INFO_updated'] = vcf_df['INFO'] + ';' +  vcf_df['INFO_add']
    vcf_df['INFO_updated'] = vcf_df['INFO_updated'].fillna(vcf_df['INFO'])   
    vcf_df = vcf_df.drop(columns=['INFO_add','INFO']).rename(columns={"INFO_updated": "INFO"})
    
    
    logging.info('Writing out the main VCF body to output file')
    # using append mode because header already written
    vcf_df.to_csv(out_file, mode='a', index=False, sep='\t', header=True)


    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse ClinVar vcf and append extra information from the ClinVar XML')
    req_grp = parser.add_argument_group(title='Required')
    req_grp.add_argument('-x', '--xml',type=str, help='ClinVar XML file (can be .gz)', required=True)
    req_grp.add_argument('-i', '--input', type=str, help="ClinVar input vcf file name (can be .gz)",required=True)
    req_grp.add_argument('-o', '--out', type=str, help="Output vcf file name (non gz)",required=True)
    parser.add_argument('-l', '--log',type=str, help="Log file name")
    
    
    args = parser.parse_args()
    
    if args.log:    
        logging.basicConfig(filename=args.log, level=logging.INFO)
    else:
        out_path = os.path.dirname(args.out)
        logging.basicConfig(filename=(out_path + "/match_alleleid.log"), level=logging.WARNING)
        
        
    log = logging.getLogger(__name__)
    
    logging.info(f"Start time: {datetime.now()}\n")
    logging.info(args)
    
    expand_clinvar_vcf(xml_file=args.xml, out_file=args.out, 
                       vcf_file=args.input)

    logging.info(f"\nEnd time: {datetime.now()}")

    

