Assembly and Mapping 
---------------------

Steps:
~~~~~~~

DAG
~~~

.. image:: Assembly_Indexing_mapping.PNG   (Try saving the DAG png with the same name as the .rst file)
   :alt: Assembly and mapping DAG

Requires
~~~~~~~~

Programs required
~~~~~~~~~~~~~~~~~~

Example of Sample File
~~~~~~~~~~~~~~~~~~~~~~

::

    Title	Paired_end_project

    #SampleID	Type	Path    lane
    Sample1	Forward	/path/to/Sample1_F1.fastq.gz 1
    Sample1	Forward	/path/to/Sample1_F2.fastq.gz 2
    Sample1	Reverse	/path/to/Sample1_R1.fastq.gz 1
    Sample1	Reverse	/path/to/Sample1_R2.fastq.gz 2
    Sample2	Forward	/path/to/Sample2_F1.fastq.gz 1
    Sample2	Reverse	/path/to/Sample2_R1.fastq.gz 1
    Sample2	Forward	/path/to/Sample2_F2.fastq.gz 2
    Sample2	Reverse	/path/to/Sample2_R2.fastq.gz 2
