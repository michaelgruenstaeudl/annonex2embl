CHANGELOG
---------
#### Version 1.0.1 (2021.03.28)
* Minor bugfix
* Updating dependency requirements
#### Version 1.0.0 (2020.03.08)
* Minor improvements to the example files
* Minor code improvements
#### Version 0.9.5 (2020.03.06)
* Added functions that check each time an interaction with a third-party server is conducted if that server is up and running; the user is notified of this check and a possible down-time of the server via messages to the screen.
* Added optional commandline parameter "qualifiername" that allows users to specifiy the name of the qualifier that contains the product information (or other information) of a sequence feature in the resulting EMBL flatfile.
#### Version 0.9.0 (2020.01.10)
* Added optional commandline parameter to specify the metadata delimiter, with the default set to a comma (",")
* Automatically replaces the mol_type "DNA" with "genomic DNA" in the ID line of each EMBL record
* Implemented tests that check if any sequence name is duplicated in either the NEXUS or the metadata file.
* Implemented a test that checks if every sequence name in the NEXUS file has a corresponding entry in the metadata file.
#### Version 0.8.5 (2019.10.24)
* Minor code improvements to the argparse interface
* Minor code improvements to the tax_check behaviour
#### Version 0.8.0 (2019.10.16)
* Added Travis-CI support
* Minor code improvements to ensure cross-platform compatibility
#### Version 0.7.0 (2019.10.15)
* Added a function that drops all sequence records smaller than 10 unambiguous nucleotides
* Minor code improvements
#### Version 0.6.5 (2019.10.11)
* Unified information provision to user in case of warnings or exceptions
* Updated example datasets
#### Version 0.6.0 (2019.10.09)
* Added a function that compresses the output flatfile upon user request
* Added a function that conducts the taxonomy check primarily via the [ENA taxonomy service](https://www.ebi.ac.uk/ena/browse/taxonomy-service)
* Added function to post-process the individual sequence records via internal functions instead of calling Bash's sed
* Updated example files to reflect capability to handle forward and reverse annotations
#### Version 0.5.0 (2019.09.10)
* Make code compatible with Python 3
* Make the application of the product lookup optional
#### Version 0.4.6 (2019.09.09)
* Added function that adds "/codon_start=1" in CDS feature, if start and stop position of feature is uncertain
* Automatically remove gene and exon annotation features than are less than 15 nt long (as they are not permissible on ENA)
* Added a function that automatically removes all sequences that consists only of Ns (or "?"s)
* Added a function that compensate for the contraction of an annotation due to the identification of an internal stop codon
* Added a function that reads in if a charset is forward or reverse
* Added a function that automatically generates a [manifest file](https://ena-docs.readthedocs.io/en/latest/cli_01.html#manifest-file-types)
* Added a function that removes the accession number from the AC line and from the ID line upon selection by user
#### Version 0.4.5 (2018.05.22)
* Added function that converts missing sections of a sequence that are longer than 2 nucleotides into a "gap"-feature
#### Version 0.4.4 (2018.03.29)
* Separate the annonex2embl from the embl2checklist function.
#### Version 0.4.3 (2018.03.23)
* Improved formatting of Python code
* Checking if sequence names in NEX-file identical to sequence ids in csv-file
* Discarding charset_ranges that are empty
* Source feature 'organelle' implemented
#### Version 0.4.2 (2017.02.01)
* All qualifier values are formatted to consist of ASCII characters only.
* If a coding region is among the sequence features, the qualifier /\trans_table/ is added to the source feature.
* Implementation of customization of DE line.
#### Version 0.3 (2017.01.25)
* Rudimentary checklist output has been implemented.
* Taxon names are compared to [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) and converted to informal names if not listed there.
#### Version 0.2 (2017.01.22)
* Leading and trailing ambiguities per sequence are removed while maintaining correct annotations.
* Modifiers without content are no longer included in the output.
* Specification of taxonomic division and sequence version has been implemented.
* A submission mode that masks the ID line and the AC line has been implemented.
