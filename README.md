# Additional NeatSeq-Flow modules

This repository includes modules that extend the functionality of NeatSeq-Flow to additional programs.

Over time, additional modules will be added to the repo to extend NeatSeq-Flow to new and existing programs not yet included.

## Installing and executing NeatSeq-Flow

Please follow the instructions in [the main NeatSeq-Flow repo](https://github.com/bioinfo-core-BGU/neatseq-flow).

## Using the additional modules repo

1. Download and unzip this repo to a location of your choice.
2. Create a [pipeline parameter file](http://neatseq-flow.readthedocs.io/en/latest/02.build_WF.html#parameter-file-definition) using the directions [specified here](http://neatseq-flow.readthedocs.io/en/latest/02.build_WF.html#parameter-file-definition).
3. Add the path to the `module_path` [global parameter](http://neatseq-flow.readthedocs.io/en/latest/02.build_WF.html#global-parameters). (If you have more than one path, you can pass several to the `module_path` as a YAML list).
4. Execute NeatSeq-Flow as [described here](http://neatseq-flow.readthedocs.io/en/latest/02.build_WF.html#execution). Don't forget to pass your pipeline parameter file with the `-p` option. 

## Adding your own modules to NeatSeq-Flow

You can easily [write your own modules](http://neatseq-flow.readthedocs.io/en/latest/06.addnew_module.html#for-the-programmer-adding-modules), to include currently unsupported programs in NeatSeq-Flow. Adding these modules is done by adding the local path to the `module_path` line in the parameter file.

Alternatively, you can add your program without creating a module, [using the `Generic` module](http://neatseq-flow.readthedocs.io/en/latest/modules/generic.html#module-step_classes.Generic.Generic). This enables quickly adding steps to your pipeline without needing to create modules. However, you will need to define the input and output locations in the parameter file. See the [documentation](http://neatseq-flow.readthedocs.io/en/latest/modules/generic.html#module-step_classes.Generic.Generic) for the generic module.  
