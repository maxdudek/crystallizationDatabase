# HWI PDB Crystallization Database

This Python package contains tools for building and searching a database of crystallization details and protein sequences from the PDB. Features include:
* Ability to pull information from the PDB and put it in a Python object
* Automatic parsing of the crystallization details into a list of compounds and concentrations
* Standardization of chemical names by a dictionary of compounds
* Application to dynamically build and improve the dictionary

The crystallization database is stored as a list of Python objects serialized to a binary file (with a `.pkl` extension*). The Python objects come from the Structure class, and contain information including PDB ID, sequences, and crystallization details. They are created by pulling the information of every structure from the PDB which contains relevant crystallization details (~134k structures as of December 2019). A “detail parsing function” is used to turn the raw plain English details in the PDB into a consistently formatted list of chemical compounds and their concentration, which is also stored as part of the Structure objects. A compound dictionary was also built to map the most commonly extracted compound names (about 680 of them) to a set of unique standardized names (~310) for consistency.  Compound names which appear more than 30 times in the PDB are recognized, and they account for 88% of all compound names which are extracted by the detail parsing function.
 
Of the ~134k structures with chemical details, over 99k are “sensible” – that is, all of the compounds which were extracted from the crystallization details are recognized by the compound dictionary. The “sensible structures” database is the database which should be used for analysis, as complete recognition by the dictionary is a sign that the detail parsing function performed as expected, and correctly extracted the reservoir compounds. The structures which are not “sensible” may be so for a few reasons:
1. They contain uncommon compounds which are not yet recognized by the dictionary – these compounds, many of which are the name of the protein itself, must be added manually, which is time consuming
2. The detail parsing function can not properly extract compounds from the details due to unexpected syntax, spacing, or punctuation
3. The details field contains a mistake, typo, or a lack of sufficient details.

For more details and results from the 2018 summer project, see the project poster in `project_poster.pdf`.

For more information on how to use the package see the [Wiki](https://github.com/maxdudek/crystallizationDatabase/wiki).

Contact Sarah Bowman at sbowman@hwi.buffalo.edu or Max Dudek at max.dudek@pitt.edu if you have any questions.

\* The database file is too large to include directly in the repository, but can be downloaded from the data repository [Zenodo](https://doi.org/10.5281/zenodo.3931013) doi: 10.5281/zenodo.3931013. Previously these data were available from the Structure directory using the github large file storage extension, but bandwidth limitations limit the number of downloads that can occur, so we have mirrored the data to the Zenodo repository.
