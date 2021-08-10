import importlib.resources as pkg_resources

from unittest import TestCase
import xml.etree.ElementTree as ET

from clinvar_vcf.clinvar_vcf_parser import *
import tests.resources


class TestModuleFunctions(TestCase):
    """
    Test the top-level module functions in clinvar_vcf_parser
    """

    def test_get_accession(self):
        """
        get_accession behaves as expected
        """
        with pkg_resources.path(tests.resources,
            'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_accession(cva)
            expected = ['SCV000020619']
            self.assertEqual(returned, expected) 

    def test_get_clinsig_status_ordered(self):
        """
        get_clinsig behaves as expected
        """
        with pkg_resources.path(tests.resources,
            'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_clinsig(cva, field='status_ordered')
            expected = ['no assertion criteria provided']
            self.assertEqual(returned, expected)

    def test_get_clinsig_last_eval(self):
        """
        get_clinsig behaves as expected
        """
        with pkg_resources.path(tests.resources,
            'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_clinsig(cva, field='last_eval')
            expected = ['2007-03-01']
            self.assertEqual(returned, expected)

    def test_get_clinsig_description(self):
        """
        get_clinsig behaves as expected
        """
        with pkg_resources.path(tests.resources,
            'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_clinsig(cva, field='description')
            expected = ['Pathogenic']
            self.assertEqual(returned, expected)

    def test_get_clinsig_comment(self):
        """
        get_clinsig behaves as expected
        """
        with pkg_resources.path(tests.resources,
            'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_clinsig(cva, field='comment')
            expected = ['']
            self.assertEqual(returned, expected)

    def test_get_origin(self):
        """
        get_origin behaves as expected
        """
        with pkg_resources.path(tests.resources,
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
        with pkg_resources.path(tests.resources,
            'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            cva = cvs.findall('./ClinVarAssertion')
            returned = get_submitdate(cva)
            expected = ['2015-08-20']
            self.assertEqual(returned, expected)

    def test_get_submitters(self):
        """
        get_submitters behaves as expected
        """
        with pkg_resources.path(tests.resources,
            'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            returned = get_submitters(cvs)
            expected = ['OMIM']
            self.assertEqual(returned, expected)

    def test_get_traits(self):
        """
        get_traits behaves as expected
        """
        with pkg_resources.path(tests.resources,
            'ClinVarFullRelease_2021-03_1_record.xml') as fpath:
            cvs = ET.parse(fpath).getroot().find('./ClinVarSet')
            returned = get_traits(cvs, field='traits')
            expected = ['BLOOD GROUP--LUTHERAN NULL'] 
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
                                'clinvar_20210302_10_records.vcf') as fpath:
            header, vcf = read_vcf(str(fpath))
            self.assertEqual(header[0], '##fileformat=VCFv4.1') 
            self.assertEqual(len(header), 27)
            self.assertEqual(vcf['#CHROM'][0], '1')
            self.assertEqual(len(vcf), 10)

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
                                'clinvar_20170404_10_records.vcf') as fpath:
            header, vcf = read_vcf(str(fpath))
            split_vcf = split_multi_vcf(vcf)
            returned = split_vcf_info(split_vcf)
            expected = {
                'RCV000162196': 0,
                'RCV000148989': 1,
                'RCV000148988': 2,
                'RCV000424799': 3,
                'RCV000422793': 4,
                'RCV000116272': 5,
                '': 6,
                'RCV000193277': 7,
                'RCV000250556': 8,
                'RCV000235037': 9,
                'RCV000116258': 10
            }
            self.assertEqual(returned, expected)
