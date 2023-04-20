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
It is required to have Python 3.9 or higher in order to run the script. The script itself was build using Python 3.10, but can sustain lower versions. In order to properly process the data and recreate the PAF file, Minimap2 is required. ```git clone https://github.com/lh3/minimap2``` can be used to download the repository. In order to be able to use the Minimap2 tool use the following command to get into the directory and creating the correct environment: ```cd minimap2 && make```. And in order to create the PAF file for the Nanopore data use ```./minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf```. Note that the reads.fa are the same file twice, there is no secondary file required here. A complete detailed guide about Minimap2 can be found [here](https://github.com/lh3/minimap2).

### Packages
The packages used for the assembler can be found within the requirements.txt file, a full description of their purpose and use is listed below.

|Name                                   | Version           |
|---                                    |-------------------|
|__future__                             | 3.10              |
|Collections                            | 3.10              |
|Os                                     | 3.10              |
|Math                                   | 3.10              |
|Subprocess                             | 3.10              |
|Argparse                               | 3.10              |
|Re                                     | 3.10              |
|Numpy                                  | 1.24.2            |
|Bio                                    | 1.81              |
|Matplotlib                             | 3.7.1             |
|Pandas                                 | 2.0.0             |

As can been seen in the table above, eleven packages were used within the script. most packages used are in version 3.10 which is the Python version that was used, this version depends solely on the version of Python installed on the system.

* __future__: This package was used within the PAF file reader, this script reforms the PAF file into a CSV file which can be loaded in by Pandas.

* Collections: This package was used within the PAF file reader, in order to create the tuples needed for the header names of the CSV file.

* Os: The Os package was used in the building of the data frame in order to find the correct path to the required information.

* Math: The Math package, just like the Os package was used within building the data frame and normalizing the data within that data frame.

* Subprocess: This package was used in order to run Minimap2 to create the PAF file so that it can later be used for the creation of the CSV file.

* Argparse: ***MAY CHANGE LATER*** This package was used for the parsing of arguments so that the user can chose what actions the program takes and which actions are skipped

* The Re regular expressions were used within the script in order to separate reads from read names and collect these efficiently.

* Numpy: The Numpy package was used for it statistical capabilities, in order to calculate things such as GC percentage.

* Bio: The Bio package can work with a large variety of data, and was used here in order to calculate multiple statistics, together with Numpy

* Matplotlib: This plotting module was used in order to make multiple visualizations such as the GC percentage of the file.

* Pandas: The Pandas module was used for the gathering of the data within the file. Most of the assembly works with Pandas.

### Scripts
There are two main scripts that run the assembler, first there is the PAF_reader.py script, this script, this script makes sure the PAF file is set to a CSV file so that the assembler can use it. Secondly the assembler.py is the final and biggest script. Here the reads are assembled into contigs, that can be analysed and used. Here the user also gets the options to visualize data.

## Usage
The script itself can be called by a simple Python statement within the command line, as long as all the correct packages are installed. Make sure that the Python version used, is above Python 3.8

### Examples
```Python Assembly.py -f [file] -c [number] -v -SR ```
With this statement the Script Assembly.py is called upon, this script runs the PAF_reader as well, since the file created there is required here. The first argument here is the file, this will be the MinION FastQ file that will need to be inputted. Secondly, the '-C' argument stands for the amount of cores the user wants to run the program on, more cores means a faster process. The last two arguments are the "-v" and the "-SR" arguments, the visualize argument can be chosen by the user to show graphs such as the GC percentage. And the "-SR" stands for "Show reads", here the user will be shown how many reads were done, the overall quality and the longest and shortest read. The only required parameter here is the "-f" argument.

|Parameter| Required | Default     | Description |
|---      |----------|-------------|-------------|
|-f       | Yes      | Is required | The MinION FastQ file that needs to be assembled. |
|-c       | No       | 4           | The Amount of cores required. |
|-v       | No       | No          | Does the user want a visualization? |
|-SR      | No       | No          | Does the user want read information? |

As can be seen in the table above only the file is required, the rest of the arguments are automatically "No" to speed up the process. Cores is set to a standard of four.

## Contact
For any questions or issues please contact the creator as specified below.
* Lisan Eisinga
  * j.l.a.eisinga@st.hanze.nl 
