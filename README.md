# Binana_Log_File_Parser

It is python-implemented command line utility for getting interatomic interactions from BINANA output log.txt file.

## Notes

* This program is distributed in the hope that it will be useful, but without any warranty.
* Python 3.6+ is required.
* Read the Input file paragraph carefully, since this program analyzes only the main types of ligand-receptor interactions. 

### Input file

As input, Binana_log_parser accepts BINANA log.txt file that contains a description
of the ligand-receptor interactions. 
    
BUT the program will work correctly if

    1)  both the LIGAND and the RECEPTOR are represented by only ONE CHAIN 
    (if not, Warning will be raised but program will be executed) 

    2) the LIGAND is SIMPLE (has one name) or represented only by amino acids
    (if not, Warning will be raised but program will be executed) 

    3) the LIGAND and the RECEPTOR have CHAIN NAMEs and these names are different
    in case both of them are represented by amino acids
    (if not, Error will be raised and program will terminate) 

### Running

The simplest way to run Binana_log_parser (when your input log.txt and binana_log_parser.py in the same directory):

    $ python3 binana_log_parser.py
    
As a result, you will see output.txt and binana_parser.log (only needed to determine where the error occurred). 

Recommended method of use:
   
    $ python3 binana_log_parser.py --input_file /path/to/input/log/file/log.txt --out_txt /path/to/out/txt/file/out.txt


### Parameters
    $ -h | --h
    
Shows help message.

    $ --input_file PATH_TO_INPUT_FILE

Gets source log file path. Default value: log.txt.

    $ --out_txt PATH_TO_OUTPUT_TXT_FILE

Gets result output txt file path. Default value: output.txt.

    $ --out_csv PATH_TO_OUTPUT_CSV_FILE

Gets result output csv file path.