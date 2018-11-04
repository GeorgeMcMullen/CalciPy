# Sample-1 Data

The provided set of data in this directory consists of the following:

- Sample-1.xlsx: A sample input file in the same format that Nikon NIS-Elements exports it
- Sample-1-Input_Processed.csv: Output of the script with data in CSV format
- Sample-1-Input_Sheet 1.pdf: Graphical representation of the output waveforms with marks for exponential decay and Tau

## NOTE: The attached sample set of data was generated to show the operation of the script and should not be used for actual scientific purposes.

The file was run using the following command line parameters. Adjust your paths as necessary.

  python2 calci.py --decaystart .85 --decayend .10 --bounds Samples/Sample-1/Sample-1-Input.xlsx 

The output PDF will show that the script can calculate appropriate results, even though the waveforms are somewhat noisy.
