ClinVar vcf parser that extracts information from ClinVar xml file.

**Software Requirements**

 * python3.6+
 * pandas
 * tox (for tests only)


**Install**

    git clone https://github.com/okon/clinvar_vcf
    cd clinvar_vcf
    pip3 install --user .

**Test**

    # (optional)
    tox

**Usage**

    clinvar_vcf_parser [-h] -x XML -i INPUT -o OUT [-l LOG] [--pre-may-2017]
