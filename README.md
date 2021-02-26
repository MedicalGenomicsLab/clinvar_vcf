ClinVar vcf parser that extracts information from ClinVar xml file.

Two versions:
 * clinvar_old_format_vcf_parser.py for old ClinVar vcf format pre May 2017. Uses RCV IDs from INFO vcf column.
 * clinvar_vcf_parser.py for current ClinVar vcf format after May 2017. Uses vcf ID column - ClinVar Variation IDs.

**Software Requirements**

 * python3.6+ (tested with 3.6.1)
 * pandas

**Usage**

    python clinvar_vcf_parser.py -x [XML_file] -i [VCF_file] -o [VCF_parsed_file] -l [log_file]
