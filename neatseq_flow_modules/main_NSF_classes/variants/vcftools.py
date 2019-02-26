# -*- coding: UTF-8 -*-
""" 
``vcftools`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running vcftools:

Can take a VCF, gunzipped VCF or BCF file as input.

Produces an output file, as specified by the output options arguments.


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* Input files in one of the following slots (for project scope):

    * ``sample_data["VCF" | "gzVCF" | "BCF"]``
    
* Input files in one of the following slots (for sample scope):

    * ``sample_data[<sample>]["variants"]["VCF" | "gzVCF" | "BCF"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* Puts output files in the following slots (for project scope): 
    ``self.sample_data["project_data"][<output type>]``

* Puts output files in the following slots (for sample scope): 
    ``self.sample_data[<sample>]["variants"][<output type>]``

.. Note:: 

    Output type is set by redirecting the required type, *i.e.* any number of the following list of types.
    
    For extracting several INFO fields, set ``--get-INFO`` to a list of INFO elements to extract (instead of passing ``--get-INFO`` several times). See examples below.
    
    If several output types are passed, each type will be created in parallel with a different vcftools script.
    
    See the vcftools manual for details (https://vcftools.github.io/man_latest.html).
    
    ``"--freq"``, ``"--freq2"``, ``"--counts"``, ``"--counts2"``, ``"--depth"``, ``"--site-depth"``, ``"--site-mean-depth"``, ``"--geno-depth"``, ``"--hap-r2"``, ``"--geno-r2"``, ``"--geno-chisq"``, ``"--hap-r2-positions"``, ``"--geno-r2-positions"``, ``"--interchrom-hap-r2"``, ``"--interchrom-geno-r2"``, ``"--TsTv"``, ``"--TsTv-summary"``, ``"--TsTv-by-count"``, ``"--TsTv-by-qual"``, ``"--FILTER-summary"``, ``"--site-pi"``, ``"--window-pi"``, ``"--weir-fst-pop"``, ``"--het"``, ``"--hardy"``, ``"--TajimaD"``, ``"--indv-freq-burden"``, ``"--LROH"``, ``"--relatedness"``, ``"--relatedness2"``, ``"--site-quality"``, ``"--missing-indv"``, ``"--missing-site"``, ``"--SNPdensity"``, ``"--kept-sites"``, ``"--removed-sites"``, ``"--singletons"``, ``"--hist-indel-len"``, ``"--hapcount"``, ``"--mendel"``, ``"--extract-FORMAT-info"``, ``"--get-INFO"``, ``"--recode"``, ``"--recode-bcf"``, ``"--12"``, ``"--IMPUTE"``, ``"--ldhat"``, ``"--ldhat-geno"``, ``"--BEAGLE-GL"``, ``"--BEAGLE-PL"``, ``"--plink"``, ``"--plink-tped"``.

.. Warning::
    
    At the moment, you can't pass more than one ``extract-FORMAT-info`` option at once. For more than one ``extract-FORMAT-info``, create more than one instance of ``vcftools``.


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "Indicates whether to use a project or sample bowtie1 index."
    "input", "vcf | bcf | gzvcf", "Type of input to use. Default: vcf"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    vcftools1:
        module: vcftools
        base: freebayes1
        script_path: /path/to/vcftools
        scope: project
        input: vcf
        redirects:
            --recode:
            --extract-FORMAT-info: GT
            --get-INFO:
                - NS
                - DB

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T. and McVean, G., 2011. **The variant call format and VCFtools**. *Bioinformatics*, 27(15), pp.2156-2158.

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_vcftools(Step):

   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        if "scope" not in self.params:
            raise AssertionExcept("You must supply a 'scope' parameter!")
        elif self.params["scope"] not in ["project","sample"]:
            raise AssertionExcept("'scope' must be 'sample' or 'project' (case sensitive)!")
                    
        if "input" not in self.params:
            self.write_warning("No 'input' parameter passed. Looking for VCF")
            self.params["input"] = "vcf"
        elif self.params["input"] not in ["vcf","bcf","gzvcf"]:
            raise AssertionExcept("'input' must be 'vcf', 'bcf' or 'gzvcf' (case sensitive)!")
        # Check that 'redirects' contains only one of the following options:
        if "redir_params" not in self.params:
            raise AssertionExcept("You must supply an output type in 'redirects'")
            
        # Checking the length of the intersection of redirects params and the possible list of output types:
        # If 0 - none passed. Use 'recode' as default
        # If >1 - error. two passed.
        # If 1 - OK.
        # Getting possible output types from get_types_index() (see definition below)
        vcftools_output_types = ["--%s"%otype for otype in list(self.get_types_index().keys())]
        # Adding '--extract-FORMAT-info' which get's treated individually
        vcftools_output_types.append("--extract-FORMAT-info")
        self.output_types = list(set(self.params["redir_params"]) & set(vcftools_output_types))

        if len(self.output_types) == 0:
            self.output_types = {"--recode" : None}
            self.write_warning("No output type passed. Using '--recode' by default")
            
        # Converting output types to external list:
        self.output_types = {key:self.params["redir_params"][key] for key in self.output_types}  # = None
        # Remove output_types from redir_params:
        self.params["redir_params"] = {key:self.params["redir_params"][key] for key in self.params["redir_params"] if key not in list(self.output_types.keys())}

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:
            if self.params["input"] == "vcf":
                try:
                    self.sample_data[sample]["vcf"]
                except KeyError:
                    raise AssertionExcept("No VCF variants file.", sample)
            elif self.params["input"] == "bcf":
                try:
                    self.sample_data[sample]["bcf"]
                except KeyError:
                    raise AssertionExcept("No BCF variants file.", sample)
            else:
                try:
                    self.sample_data[sample]["gzVCF"]
                except KeyError:
                    raise AssertionExcept("No 'gzVCF' (gzipped VCF) variants file.", sample)

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass
        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        types_index = self.get_types_index()

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:      # Getting list of samples out of samples_hash

            sample_title = sample if sample != "project_data" else self.sample_data["Title"]

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            if self.params["input"] == "vcf":
                input_file = self.sample_data[sample]["vcf"]
            elif self.params["input"] == "bcf":
                input_file = self.sample_data[sample]["bcf"]
            else:
                input_file = self.sample_data[sample]["gzVCF"]

            output_prefix = "{dir}{step_name}_{sample}".format(dir = use_dir,
                                                               step_name = self.get_step_name(),
                                                               sample = sample_title)

            for output_type in self.output_types:
                # Name of specific script:
                self.spec_script_name = self.jid_name_sep.join([self.step,
                                                  self.name,
                                                  output_type.lstrip("-"),
                                                  sample_title])
                self.script = ""

                # Add current output_type to redir_params (letting neatseq flow deal with duplicating
                # lists and removing 'None's.)
                self.params["redir_params"][output_type] = self.output_types[output_type]
                # Get constant part of script:
                self.script += self.get_script_const()
                # Reference file:
                self.script += "--{flag} {file} \\\n\t".format(flag = self.params["input"], file = input_file)
                self.script += "--out %s \\\n\t" % output_prefix

                # # Remove the output_type from redir_params:
                del self.params["redir_params"][output_type]

                # TODO: CHECK
                for type in types_index:
                    if output_type == "--%s" % type:
                        for outype in types_index[type]:
                            self.sample_data[sample][outype] = "{prefix}.{type}.{outype}".\
                                                                    format(prefix=output_prefix,
                                                                            type=type,
                                                                            outype=outype)
                            self.stamp_file(self.sample_data[sample][outype])

                if output_type == "--extract-FORMAT-info":
                    self.sample_data[sample]["extract-FORMAT-info"] = \
                        "{prefix}.{info}.{suffix}".format(prefix = output_prefix, \
                                                          info = self.output_types["--extract-FORMAT-info"],\
                                                          suffix = "FORMAT")
                    self.stamp_file(self.sample_data[sample]["extract-FORMAT-info"])

                self.local_finish(use_dir,sample_dir)
                self.create_low_level_script()

    def get_types_index(self):
        """

        :return: Dict connecting output type with output file suffixes.
        """
        return {"recode":["vcf"],
                "recode-bcf":["bcf"],
                "012":["012","012_indv","012_pos"],
                "IMPUTE":["impute.hap","impute.hap.legend","impute.hap.indv"],
                "ldhat":["ldhat.sites","ldhat.locs"],
                "ldhat-geno":["ldhat.sites","ldhat.locs"],
                "BEAGLE-GL":["BEAGLE.GL"],
                "BEAGLE-PL":["BEAGLE.PL"],
                "plink":["ped","map"],
                "plink-tped":["tped","tfam"],
                "freq":["freq"],
                "freq2":["freq2"],
                "counts":["counts"],
                "counts2":["counts2"],
                "depth":["depth"],
                "site-depth":["site-depth"],
                "site-mean-depth":["site-mean-depth"],
                "geno-depth":["geno-depth"],
                "hap-r2":["hap-r2"],
                "geno-r2":["geno-r2"],
                "geno-chisq":["geno-chisq"],
                "hap-r2-positions":["hap-r2-positions"],
                "geno-r2-positions":["geno-r2-positions"],
                "interchrom-hap-r2":["interchrom-hap-r2"],
                "interchrom-geno-r2":["interchrom-geno-r2"],
                "TsTv":["TsTv"],
                "TsTv-summary":["TsTv-summary"],
                "TsTv-by-count":["TsTv-by-count"],
                "TsTv-by-qual":["TsTv-by-qual"],
                "FILTER-summary":["FILTER-summary"],
                "site-pi":["site-pi"],
                "window-pi":["window-pi"],
                "weir-fst-pop":["weir-fst-pop"],
                "het":["het"],
                "hardy":["hardy"],
                "TajimaD":["TajimaD"],
                "indv-freq-burden":["indv-freq-burden"],
                "LROH":["LROH"],
                "relatedness":["relatedness"],
                "relatedness2":["relatedness2"],
                "site-quality":["site-quality"],
                "missing-indv":["missing-indv"],
                "missing-site":["missing-site"],
                "SNPdensity":["SNPdensity"],
                "kept-sites":["kept-sites"],
                "removed-sites":["removed-sites"],
                "singletons":["singletons"],
                "hist-indel-len":["hist-indel-len"],
                "hapcount":["hapcount"],
                "mendel":["mendel"],
                "get-INFO":["get-INFO"]
                }
