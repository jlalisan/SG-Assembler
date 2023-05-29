# An inside look into a genome assembler #
The Study Genome assembler (SG-assembler) is a versatile assembler, that can assemble the smallest files into contigs for the genome and visualize the data within the PacBio MinION file. Because the assembler was fully written in Python 3.10 and documented it is great or teaching students about genome assembly and the inner workings of an assembler.

## Table of content
- [Installation](#installation)
    * [Prerequisites](#prerequisites)
    * [Packages](#packages)
    * [Scripts](#scripts)
- [Usage](#usage)
    * [Examples](#examples)
- [Contact](#contact)

## Installation
In order to run the script on a data set of choice the repository needs to be cloned. This can be done by either downloading the repository from GitHub or by entering the following command into a command line. ```git clone https://github.com/jlalisan/SG-Assembler``` After the repository has been cloned, the required dependencies can be found within the requirements.txt file and can be used with the following command: ```pip install -r requirements.txt``` This command needs to be used within the correct location, where the script is located.

### Prerequisites
It is required to have Python 3.8 or higher in order to run the script. The script itself was build using Python 3.10, but can sustain lower versions. In order to properly process the data and recreate the PAF file, Minimap2 is required. ```git clone https://github.com/lh3/minimap2``` can be used to download the repository. In order to be able to use the Minimap2 tool use the following command to get into the directory and creating the correct environment: ```cd minimap2 && make```. And in order to create the PAF file for the Nanopore data use ```./minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf```. Note that the reads.fa are the same file twice, there is no secondary file required here (own file required). A complete detailed guide about Minimap2 can be found [here](https://github.com/lh3/minimap2).

### Packages
The packages used for the assembler can be found within the requirements.txt file, a full description of their purpose and use is listed below.

|Name                                   | Version           |
|---                                    |-------------------|
|BioPython                              | 1.81              |
|Docopt                                 | 0.6.2             |
|Logging                                | 3.10              |
|Networkx                               | 3.1               |
|Os                                     | 3.10              |
|Pandas                                 | 2.0.1             |
|Subprocess                             | 3.10              |
|Sys                                    | 3.10              |
|Time                                   | 3.10              |


As can been seen in the table above, nine packages were used within the script. most packages used are in version 3.10 which is the Python version that was used, this version depends solely on the version of Python installed on the system.

Within the script the packages were used for:

* BioPython: Used for the reading of the FastQ file and the ordering of the sequences.

* Docopt, used for the building of a command line interface.

* Logging, used for the printing of statements.

* Networkx, used for the building of the multidigraph and the removal of isolated nodes.

* Os, used to check if paths exist for files.

* Pandas, used to create a dataframe from the PAF file.

* Subprocess, used to use minimap2 in case a PAF file does not exist.

* Sys, used for the exiting of programs and setting of recursion limit

* Time, used to calculate how long the assembler takes to run.


### Scripts
There are two main scripts that run the assembler, first there is the Assembly script, This script runs the entire assembly, however if the user choses too the script makes use of 'visualise.py' which give visualisations of the FastQ file.
## Usage
The script itself can be called by a simple Python statement within the command line, as long as all the correct packages are installed. Make sure that the Python version used, is above Python 3.8

### Examples
```python Assembly.py <fastq_file> [--paf_file=<paf_file>] [--output_file=<output_file>] [--graphics]```

With this statement the Script Assembly.py is called upon, this script runs the visualization as well. The first argument here is the file, this will be the MinION FastQ file that will need to be inputted, this is also the only required file. The second argument is the paf_file argument, that will be created if it is left empty, similar to the third argument which is the output file for the assembly, if this is left empty this will also be created based on the FastQ name. Lastly the graphics argument is for the visualization of the FastQ file if called for.

|Parameter      | Required | Default     | Description        |
|---            |----------|-------------|--------------------|
|<fastq_file>   | Yes      | Is required | FastQ MinION file. |
|[--paf_file]   | No       | fastq.paf   | Overlaps file.     |
|[--graphics]   | No       | No          | Visualization?     |
|--[output_file]| No       | fastq.fasta | Output file name   |

As can be seen in the table above only the file is required, the rest of the arguments are automatically "No" to speed up the process.

## Contact
For any questions or issues please contact the creator as specified below.
* Lisan Eisinga
  * j.l.a.eisinga@st.hanze.nl 
