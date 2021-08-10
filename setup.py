from setuptools import setup

setup(
    entry_points={
        'console_scripts': 'clinvar_vcf_parser = clinvar_vcf.clinvar_vcf_parser:main'
    }
)
