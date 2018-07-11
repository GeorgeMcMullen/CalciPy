# CalciPy: A Python script for processing ratiometric calcium fluorescence data
One of the primary challenges in processing ratiometric data is that the two signals which make up the ratiometric data are actually captured separately, and it is not known which is which. CalciPy processes interpolated calcium fluorescence wavelength data, identifies the proper ratio sequence, and calculates metrics from the resulting signal. See below for more detail. While there is a variety of software for processing calcium fluorescence data, each of them have one or more shortcomings, such as:

- Few solutions actually process ratiometric data provided by two wavelengths of light
- Some require manual effort to determine the proper ratio and subsequent calculations
- Those that can process ratiometric data automatically are limited (e.g. only able to process a single cell line at a time)
- Finally, they are not widely available, are closed source, or otherwise proprietary limiting the auditability of the experiments

CalciPy's goal is to enable high throughput, yet customizable, processing of ratiometric calcium fluorescence data, in an open and autidible fashion. It can also serve as a starting point for processing data that is exported from microscopy imaging software, such as Nikon's NIS-Elements.

## Requirements, Dependencies, and Installation
CalciPy expects that the input file is an Excel file containing ratiometric calcium fluorescent data, which is in the same structure that software packages like Nikon NIS-Elements exports data. See the samples directory for more information.

The script requires Python (v2.7 at least). Furthermore, the following Python libraries are required and can be installed with Pip.

- openpyxl
- numpy                                                                      
- scipy                           
- matplotlib                      
- analytic-wfm

Other than that, the script can be run as is when downloaded as a signle file.                  

## Running / Command Line Arguments
CalciPy can be run simply by providing the path to the input file as a command line argument.

    usage: calci.py [-h] [-O OUTPUTDIR] [--sheetname SHEETNAME]
                    [--sheetnum SHEETNUM] [--column COLUMN]
                    [--bgreduce {average,point,none}] [--lookahead LOOKAHEAD]
                    [--delta DELTA] [--invert] [--limit LIMIT]
                    [--decaystart DECAYSTART] [--decayend DECAYEND] [--ymax YMAX]
                    [--ymin YMIN] [--show]
                    filename
    
    positional arguments:
      filename              input filename to process
    
    optional arguments:
      -h, --help            show this help message and exit
      -O OUTPUTDIR, --outputdir OUTPUTDIR
                            output files to a directory
      --sheetname SHEETNAME
                            process only a specific worksheet by name (use quotes)
      --sheetnum SHEETNUM   process only a specific worksheet by number
      --column COLUMN       process only a specific cell line column of data
                            (typically 1-10)
      --bgreduce {average,point,none}
                            method for reducing background data
      --lookahead LOOKAHEAD
                            peak detection look-ahead parameter (default=30)
      --delta DELTA         peak detection delta parameter (default=0)
      --invert              invert the waveform (special cases only)
      --limit LIMIT         limit processing to the first X wavelets (default=0 as
                            disabled)
      --decaystart DECAYSTART
                            the relative amplitude from the peak where the
                            analysis will start
      --decayend DECAYEND   the relative amplitude from the peak where the
                            analysis will end
      --ymax YMAX           Y-axis maximum
      --ymin YMIN           Y-axis minimum
      --show                show the graph (with matplotlib)

## Methods
Raw data is exported from software such as Nikon NIS-Elements into an Excel compatible file. See the samples directory for an example data set and format. As mentioned, the two signals which make up the ratiometric data are actually captured separately, represented as interpolated rows in a single dataset, and it is not defined which set of interpolated rows corresponds to what part of the ratio.

![1 - Original Data](Documentation/Images/1%20-%20Original%20Data.png "1 - Original Data")

To overcome this challenge, the script will analyze different combinations of the data, comparing it with an expected sharp rise and slow(er) decay, typical of calcium fluorescence. The script does this by first iterating over the captured cell lines and subtracting background noise (which should also be captured as the last cell line/column).

![2 - Background Signal Reduced](Documentation/Images/2%20-%20Background%20Signal%20Reduced.png "2 - Background Signal Reduced")

Next it de-interpolates the rows to capture the individual signals. At this stage it is not known which signal represents which wavelength of data.

![3a - De-interpolated Signal Data (Unknown Wavelength)](Documentation/Images/3a%20-%20De-interpolated%20Signal%20Data%20(Unknown%20Wavelength).png "3a - De-interpolated Signal Data (Unknown Wavelength)") ![3b - De-interpolated Signal Data (Unknown Wavelength)](Documentation/Images/3b%20-%20De-interpolated%20Signal%20Data%20(Unknown%20Wavelength).png "3b - De-interpolated Signal Data (Unknown Wavelength)") 

It will divide these two separate signals by eachother to calculate two ratiometric signals, which are purely the inverse of eachother.

![4a - Ratiometric Signal (Chosen)](Documentation/Images/4a%20-%20Ratiometric%20Signal%20(Chosen).png "4a - Ratiometric Signal (Chosen)") ![4b - Ratiometric Signal (Discarded)](Documentation/Images/4b%20-%20Ratiometric%20Signal%20(Discarded).png "4b - Ratiometric Signal (Discarded)") 

The script then uses look-ahead peak detection to determine local maxima and minima. These detected peaks are used to calculate rise and decay times and identify which combination of rows presents the correct ratio.

![6 - Peak Detection](Documentation/Images/6%20-%20Peak%20Detection.png "6 - Peak Detection")

A centered window moving average is applied to the decay data only (each individual set of maximum to minimum), slightly smoothing the signal to obtain more consistent initial values for curve fitting. Python’s curve fitting functions found in SciPy are then used with an exponential decay function to establish the signal’s parameters, most significantly the decay rate (which can also be used to calculate Tau).

![7 - Curve Fitting (Amplitude Range 85-10)](Documentation/Images/7%20-%20Curve%20Fitting%20(Amplitude%20Range%2085-10).png "7 - Curve Fitting (Amplitude Range 85-10)")

The original signal, peaks, fitted curve results, and marker for the position of Tau are then plotted and exported as PDFs for easy visual review.

![8 - Tau Marker Placement](Documentation/Images/8%20-%20Tau%20Marker%20Placement.png "8 - Tau Marker Placement")

In addition, calculated results for rate constant, Tau, decay time, amplitude, diastolic, beat rate, beat variation, and velocity are all outputted to a comma separated value text file (CSV) for further analysis. The script can be configured with different peak detection look-ahead windows, amplitude windows, background noise reduction, smoothing filters and more, and has shown consistent results even when run against noisy or arrhythmic signals. 

## Caveats and Troubleshooting
This software comes without any warrantee or guarantee. While it has shown to provide consistent results, even when run against noisy or arrhythmic signals, verification is still encouraged. Below are some notes of what to look out for during the verification process.

- Use the --show option or view the outputted PDF files to visually check the results
- Data that is extremely noisy may still pose problems with analysis. You should expect this if you are unable to see visible peaks in the data manually.
- If you see that the larger peaks are detected along side very small peaks, test out different --lookahead and --delta parameters
- The same should be true for if the script is not able to provide the correct ratio. To force the script to invert the ratios, use the --invert option.
- High precision microscopy instruments can sometimes lose focus or become otherwise uncalibrated over time, due to environmental factors such as temperature or humidity. Use the --limit option to only process a subset of the data.
- Some code is not implemented in the most performance or elegant manner possible. One reason for this is so that is a little easier those who are new to Python and haven't learned all its powerful capabilities with handling data.

## Special Thanks
I'd like to thank the following people at the Stanford Cardiovascular Institute who provided data, feedback, and guidance while this script was being written.

- Rajani Shrestha, Research Associate
- Timon Seeger, MD / Postdoctoral Fellow
- Lam C. Keung, PhD
- Ioannis Karakikes, PhD

