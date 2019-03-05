#from distutils.core import setup
from setuptools import find_packages, setup

setup(
    name                = 'NeatSeq-Flow-modules',
    version             = '1.6.0',
    author              = 'Menachem Sklarz',
    author_email        = 'sklarz@bgu.ac.il',
    maintainer          = 'Menachem Sklarz',
    maintainer_email    = 'sklarz@bgu.ac.il',
    url                 = 'http://neatseq-flow.readthedocs.io/projects/neatseq-flow-modules/en/latest/',
    description         = 'Shared modules and workflows for use with NeatSeq-Flow',
    license             = 'Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description    =  open('README').read(),
    download_url        = 'https://github.com/bioinfo-core-BGU/neatseq-flow-modules',
    platforms           = ["POSIX","Windows"],
    packages            = find_packages(),
    include_package_data= True,  # See  MANIFEST.in
    data_files          = [('NeatSeq-Flow-Workflows',['Workflows/Assembly_Indexing_mapping.yaml',
                                                      'Workflows/BLAST_db.yaml',
                                                      'Workflows/BLAST_fasta.yaml',
                                                      'Workflows/ChIP_seq.yaml',
                                                      'Workflows/Clustering.yaml',
                                                      'Workflows/Metagenomics.yaml',
                                                      'Workflows/parameter_set_QIIME/qiime1_mapping.txt',
                                                      'Workflows/parameter_set_QIIME/qiime1_params.txt',
                                                      'Workflows/parameter_set_QIIME/QIIME_workflow.yaml',
                                                      'Workflows/parameter_set_QIIME/sample_file_paired_small.nsfs',
                                                      'Workflows/RNA_seq_reference.yaml',
                                                      'Workflows/RNA_seq_Trinity.yaml',
                                                      'Workflows/variant_calling.yaml']),
                            ('NeatSeq-Flow-Workflows/Sample_sets',['Workflows/sample_sets/ChIP_tabular.nsfs',
                                                                   'Workflows/sample_sets/Fasta.nsfs',
                                                                   'Workflows/sample_sets/PE.nsfs',
                                                                   'Workflows/sample_sets/PE_ChIP.nsfs',
                                                                   'Workflows/sample_sets/PE_tabular.nsfs',
                                                                   'Workflows/sample_sets/SE.nsfs',])],
    # install_requires    = [],
    classifiers         = [
                          'Development Status :: 4 - Beta',
                          'Environment :: Console',
                          'Intended Audience :: End Users',
                          'Intended Audience :: Developers',
                          'License :: OSI Approved :: Python Software Foundation License',
                          'Operating System :: Microsoft :: Windows',
                          'Operating System :: POSIX',
                          'Programming Language :: Python',
                          ],
    )
    

    
