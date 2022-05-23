import filecmp
import importlib.resources as pkg_resources
import os.path
import shlex
import sys
import tempfile

from pathlib import Path
from unittest import TestCase, skipIf
import xml.etree.ElementTree as ET

from clinvar_vcf.clinvar_vcf_parser import (
    expand_clinvar_vcf,
    get_accession,
    get_clinsig,
    get_origin,
    get_submitdate,
    get_submitters,
    get_traits,
    join_entries,
    read_vcf,
    replace_sep,
    remove_newlines_and_tabs,
    split_multi_vcf,
    split_vcf_info,
    main)

import tests.resources


CV2021_XML = '/reference/data/clinvar/xml/ClinVarFullRelease_2021-03.xml.gz'
CV2021_VCF = '/reference/data/clinvar/vcf_GRCh37/clinvar_20210302.vcf.gz'
CV2021_OUT = ('/working/genomeinfo/share/olgaK/parse_ClinVar_vcfs/'
              'parsed_vcf_2021Mar/clinvar_20210302.xml_parsed.vcf')

CV2021_INTEGRATION_TEST_RESOURCES = all(
    os.path.isfile(f) for f in (CV2021_XML, CV2021_VCF, CV2021_OUT))


class TestModuleFunctions(TestCase):
    """
    Test the top-level module functions in clinvar_vcf_parser
    """

    def test_expand_clinvar_vcf(self):
        with (
            tempfile.TemporaryDirectory() as tempdir,
            pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as xmlpath,
            pkg_resources.path(
                tests.resources,
                'clinvar_20210302_1_record.vcf') as vcfpath):
            outfile = Path(tempdir) / 'parsed.vcf'
            expand_clinvar_vcf(str(xmlpath), str(vcfpath), outfile)
            returned = outfile.read_text()
            expected = pkg_resources.read_text(
                tests.resources,
                'clinvar_20210302_1_record_annotated.vcf')
            self.assertEqual(returned, expected)

    def test_get_accession(self):
        """
        get_accession behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_accession(cva)
            expected = ['SCV001214463']
            self.assertEqual(returned, expected)

    def test_get_clinsig_status_ordered(self):
        """
        get_clinsig behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_clinsig(cva, field='status_ordered')
            expected = ['criteria provided, single submitter']
            self.assertEqual(returned, expected)

    def test_get_clinsig_last_eval(self):
        """
        get_clinsig behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_clinsig(cva, field='last_eval')
            expected = ['2019-12-11']
            self.assertEqual(returned, expected)

    def test_get_clinsig_description(self):
        """
        get_clinsig behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_clinsig(cva, field='description')
            expected = ['Uncertain significance']
            self.assertEqual(returned, expected)

    def test_get_clinsig_comment(self):
        """
        get_clinsig behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_clinsig(cva, field='comment')
            expected = ['This sequence change replaces alanine with threonine at codon 36 of the SAMD11 protein (p.Ala36Thr). The alanine residue is weakly conserved and there is a small physicochemical difference between alanine and threonine. The frequency data for this variant in the population databases is considered unreliable, as metrics indicate poor data quality at this position in the ExAC database. This variant has not been reported in the literature in individuals with SAMD11-related conditions. Algorithms developed to predict the effect of missense changes on protein structure and function output the following: SIFT: "Tolerated"; PolyPhen-2: "Benign"; Align-GVGD: "Class C0". The threonine amino acid residue is found in multiple mammalian species, suggesting that this missense change does not adversely affect protein function. These predictions have not been confirmed by published functional studies and their clinical significance is uncertain. In summary, the available evidence is currently insufficient to determine the role of this variant in disease. Therefore, it has been classified as a Variant of Uncertain Significance.']
            self.assertEqual(returned, expected)

    def test_get_origin(self):
        """
        get_origin behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_origin(cva, as_set=False)
            expected = ['germline']
            self.assertEqual(returned, expected)

    def test_get_submitdate(self):
        """
        get_submitdate behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_submitdate(cva)
            expected = ['2020-02-06']
            self.assertEqual(returned, expected)

    def test_get_submitters(self):
        """
        get_submitters behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            returned = get_submitters(cvs)
            expected = ['Invitae']
            self.assertEqual(returned, expected)

    def test_get_traits(self):
        """
        get_traits behaves as expected
        """
        with pkg_resources.path(
                tests.resources,
                'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            returned = get_traits(cvs)
            expected = ['not provided']
            self.assertEqual(returned, expected)

    def test_join_entries(self):
        """
        join_entries works as expected
        """
        dct = {'a': ['1', '2', '3'], 'b': {'4', '5', '6'}}
        returned = join_entries(dct)
        expected = {'a': '1|2|3', 'b': '4|5|6'}
        self.assertEqual(returned, expected)

    def test_read_vcf(self):
        """
        read_vcf returns (header, vcf)
        """

        with pkg_resources.path(tests.resources,
                                'clinvar_20210302_1_record.vcf') as fpath:
            header, vcf = read_vcf(str(fpath))
            self.assertEqual(header[0], '##fileformat=VCFv4.1')
            self.assertEqual(len(header), 27)
            self.assertEqual(vcf['#CHROM'][0], '1')
            self.assertEqual(len(vcf), 1)

    def test_remove_newlines_and_tabs(self):
        """
        remove_newlines_and_tabs behaves as expected
        """
        string = 'foo\nbar\tbaz\rbiz'
        returned = remove_newlines_and_tabs(string)
        expected = 'foo bar baz biz'
        self.assertEqual(returned, expected)

    def test_replace_str(self):
        """
        replace_sep behaves as expected when sep and replace_with are lists
        """
        string = 'a;b|c d"e,f'
        sep = [';', '|', ' ', '"', ',']
        replace_with = [':', ':', '_', '', '&']
        returned = replace_sep(string, sep, replace_with)
        expected = 'a:b:c_de&f'
        self.assertEqual(returned, expected)

    def test_split_multi_vcf(self):
        """
        split_multi_vcf works as expected on old-format vcf
        """
        # record 7 has multiple ALT: A,C
        with pkg_resources.path(tests.resources,
                                'clinvar_20170404_10_records.vcf') as fpath:
            header, vcf = read_vcf(str(fpath))
            self.assertTrue(len(vcf) == 10)
            self.assertTrue(vcf['ALT'][6] == 'A,C')
            self.assertEqual(vcf['INFO'][6], 'RS=201073369;RSPOS=955619;dbSNPBuildID=137;SSR=0;SAO=0;VP=0x050128000a05040026000100;GENEINFO=AGRN:375790;WGT=1;VC=SNV;PM;PMC;SLO;NSM;REF;ASP;VLD;KGPhase3;CLNALLE=2;CLNHGVS=NC_000001.10:g.955619G>C;CLNSRC=.;CLNORIGIN=1;CLNSRCID=.;CLNSIG=255;CLNDSDB=MedGen;CLNDSDBID=CN169374;CLNDBN=not_specified;CLNREVSTAT=conf;CLNACC=RCV000193277.2;CAF=0.9912,.,0.008786;COMMON=1')

            split_vcf = split_multi_vcf(vcf)
            self.assertTrue(len(split_vcf) == 11)
            self.assertTrue(split_vcf['ALT'][6] == 'A')
            self.assertTrue(split_vcf['ALT'][7] == 'C')
            self.assertEqual(split_vcf['INFO'][6], 'RS=201073369;RSPOS=955619;dbSNPBuildID=137;SSR=0;SAO=0;VP=0x050128000a05040026000100;GENEINFO=AGRN:375790;WGT=1;VC=SNV;PM;PMC;SLO;NSM;REF;ASP;VLD;KGPhase3;CLNALLE=1;CLNHGVS=.;CLNSRC=.;CLNORIGIN=.;CLNSRCID=.;CLNSIG=.;CLNDSDB=.;CLNDSDBID=.;CLNDBN=.;CLNREVSTAT=.;CLNACC=.;CAF=0.9912,.,0.008786;COMMON=1')
            self.assertEqual(split_vcf['INFO'][7], 'RS=201073369;RSPOS=955619;dbSNPBuildID=137;SSR=0;SAO=0;VP=0x050128000a05040026000100;GENEINFO=AGRN:375790;WGT=1;VC=SNV;PM;PMC;SLO;NSM;REF;ASP;VLD;KGPhase3;CLNALLE=1;CLNHGVS=NC_000001.10:g.955619G>C;CLNSRC=.;CLNORIGIN=1;CLNSRCID=.;CLNSIG=255;CLNDSDB=MedGen;CLNDSDBID=CN169374;CLNDBN=not_specified;CLNREVSTAT=conf;CLNACC=RCV000193277.2;CAF=0.9912,.,0.008786;COMMON=1')

    def test_split_vcf_info(self):
        """
        split_vcf_info works as expected on old-format vcf
        """
        with pkg_resources.path(tests.resources,
                                'clinvar_20170404_13_records.vcf') as fpath:
            header, vcf = read_vcf(str(fpath))
            split_vcf = split_multi_vcf(vcf)
            returned = split_vcf_info(split_vcf)
            expected = {
                'RCV000162196': {0},
                'RCV000148989': {1},
                'RCV000148988': {2},
                'RCV000424799': {3},
                'RCV000422793': {4},
                'RCV000116272': {5},
                '': {13,6},
                'RCV000193277': {7},
                'RCV000250556': {8},
                'RCV000235037': {9},
                'RCV000116258': {10},
                'RCV000224503': {12},
                'RCV000159614': {11,14,15},
                'RCV000192572': {11,14,15},
                'RCV000286085': {11,14,15}
                
            }
            self.assertEqual(returned, expected)


class TestMain(TestCase):
    """
    Integration tests
    """

    @skipIf(not CV2021_INTEGRATION_TEST_RESOURCES,
            'Runs only when integration test resources are available')
    def test_cv2021(self):
        argv = sys.argv[:]
        try:
            with tempfile.TemporaryDirectory() as tempdir:
                outfile = Path(tempdir) / 'parsed.vcf'
                logfile = Path(tempdir) / 'parsed.vcf.log'
                sys.argv = shlex.split(
                    f'clinvar_vcf_parser -x {CV2021_XML} -i {CV2021_VCF} '
                    f'-o {outfile} -l {logfile}')
                main()
                if not filecmp.cmp(outfile, CV2021_OUT, shallow=False):
                    diffmsg = [f'{outfile} differs from expected {CV2021_OUT}:']
                    found = False
                    with open(outfile) as f1, open(CV2021_OUT) as f2:
                        for n, (line1, line2) in enumerate(zip(f1, f2), 1):
                            if line1 != line2:
                                found = True
                                diffmsg.append(f'first difference at line {n}:')
                                diffmsg.append(f' - {line2}')
                                diffmsg.append(f' + {line1}')
                        if not found:
                            f1.seek(0)
                            f2.seek(0)
                            ll1, ll2 = len(f1.readlines()), len(f2.readlines())
                            assert ll1 != ll2
                            diffmsg.append(
                                f'first {n} lines match but {outfile} has '
                                f'{ll1} lines and {CV2021_VCF} has {ll2}')
                    self.fail('\n'.join(diffmsg))
        finally:
            sys.argv = argv
