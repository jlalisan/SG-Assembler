# An inside look into a genome assembler #
Small Long-read Assembler for Educational Purposes (SLAEP) is a versatile assembler, that can assemble the smallest data sets into contigs for the genome and visualize the data within the PacBio MinION data set. Because the assembler was fully written in Python 3.10 and documented it is great or teaching students about genome assembly and the inner workings of an assembler.

## Table of content
- [Installation](#installation)
    * [Prerequisites](#prerequisites)
    * [Packages](#packages)
    * [Scripts](#scripts)
- [Usage](#usage)
    * [Examples](#examples)
- [Contact](#contact)

## Installation
In order to run the script on a data set of choice the repository needs to be cloned. This can be done by either downloading the repository from GitHub or by entering the following command into a command line. ```git clone https://github.com/jlalisan/SLAEP``` After the repository has been cloned, the required dependencies can be found within the requirements.txt file and can be used with the following command: ```pip install -r requirements.txt``` This command needs to be used within the correct location, where the script is located.

### Prerequisites
It is required to have Python 3.8 or higher in order to run the script. The script itself was build using Python 3.10, but can sustain lower versions. In order to properly process the data and recreate the PAF data set, Minimap2 is required. ```git clone https://github.com/lh3/minimap2``` can be used to download the repository. In order to be able to use the Minimap2 tool use the following command to get into the directory and creating the correct environment: ```cd minimap2 && make```. And in order to create the PAF data set for the Nanopore data use ```./minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf```. Note that the reads.fa are the same data set twice, there is no secondary data set required here. A complete detailed guide about Minimap2 can be found [here](https://github.com/lh3/minimap2).

### Packages
The packages used for the assembler can be found within the requirements.txt data set, a full description of their purpose and use is listed below.

| Python Library | Version | Usage within the script         |
| -------------- | --------| --------------------------------|
| bioPython      | 1.81    | Parsing the FastQ data set.     |
| collections    | 3.10    | Creating visualization.         |
| cProfile       | 3.10    | Checking memory efficiency.     |
| docopt         | 0.6.2   | Creating command line interface.|
| logging        | 3.10    | Printing of statements.         |
| matplotlib     | 3.7.1   | Visualizing the plots.          |
| networkx       | 3.1     | Making overlap graph.           |
| numpy          | 1.24.3  | Calculating complexity.         |
| os             | 3.10    | Finding paths to the data set.  |
| pandas         | 2.0.1   | Parsing PAF data set.           |
| psutil         | 5.9.5   | Checks run time and efficiency. |
| scipy          | 1.10.1  | Plotting functions.             |
| seaborn        | 0.12.1  | Improving plot quality.         |
| subprocess     | 3.10    | Calling on Minimap2.            |
| sys            | 3.10    | Setting recursion limit.        |
| unittest       | 3.10    | Building tests for script.      |

Seventeen packages were used for the script, each package was used for a specific reason; 

* BioPython: Package was used for the parsing of the FastQ data set, and to put the reads in a dictionary with their read names as an ID.

* Collections: This package is used to collect all the bases from the FastQ data set for visualization.

* cProfile: Used for the checking of the performance and run time of the script. Outputs a performance.prof data set.

* Docopt: This package was used for the parsing of arguments so that the user can chose what actions the program takes and which actions are skipped.

* Logging: Gives the user the messages from the script clearly legible so they wont be missed.

* Matplotlib: This plotting module was used in order to make multiple visualizations such as the GC percentage of the data set.

* NetworkX: Creates the multi directed graph for traversal filled with the reads from the FastQ as nodes and the PAF data as edges.

* Numpy: The Numpy package was used for it statistical capabilities, in order to calculate things such as GC percentage.

* Os: The Os package was used in the building of the data frame in order to find the correct path to the required information.

* Pandas: The Pandas module was used for the gathering of the data within the PAF data set.

* Psutil: Checks the memory usage of the script and outputs this to the user.

* Scipy: Works with seaborn and matplotlib in order to enhance the quality of the plots for the FastQ data set. 

* Seaborn: Gives the plots a proper background and enhances figure quality.

* Subprocess: This package was used in order to run Minimap2 to create the PAF data set.

* Sys: Used to check if the recursion limit can be set, and sets it if needed and possible.

* Unittest: Used to create a unit test for all functions so the output can be checked against expected output in case the code is being changed for educational purposes.


### Scripts
There are three scripts in the assembler, there is the 'assembly.py' which runs the entire assembly, with the option '--graphics' the user can call upon the usage of the visualization script which visualizes aspects of the FastQ data set in order to see if an assembly would be accurate. Lastly there is the "unittest.py", this script is available in the case of code changes for educational purposes and the tests can be adjusted to fit specific needs.

## Usage
The script itself can be called by a simple Python statement within the command line, as long as all the correct packages are installed. Make sure that the Python version used, is above Python 3.8

### Examples
```Python Assembly.py <fastq_file> [--paf_file=<paf_file>] [--outputS_file=<output_file>] [--graphics]```
This statement activates the 'assembly.py' script, this runs the visualization as well if requested. The first argument is the FastQ data set, this is a required argument, and requires a MinION data set, an example data set 'foo.fq' is present in the repository for testing. The second and third argument are not required. the PAF data set is created with MiniMap2 if it is not presented with the parameters, and the output data set is created from the FastQ name if not specified. The graphics can be shown with the visualization script only if the option is selected.

| Parameter   | Required? | Description                        | Default |
| ----------- | --------- | ---------------------------------- | ------- |
| fastq_file  | Yes       | The FastQ data set for assembly    | Required|
| paf_file    | No        | PAF data set of the FastQ data set | None    | 
| output_file | No        | Output data set for the assembly   | None    |
| graphics    | No        | Visualize the FastQ data set       | None    |

The parameters are listed above for usage. A complete run with the test data set 'foo.fq' can be done, this data however contains dummy data, which means the visualizations will not be optimal since the quality score for each base has been set to 'I'.

## Contact
For any questions or issues please contact the creator as specified below.
* Lisan Eisinga
  * j.l.a.eisinga@st.hanze.nl 
