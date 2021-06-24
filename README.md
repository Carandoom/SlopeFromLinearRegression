# Slope From Linear Regression [![DOI](https://zenodo.org/badge/373492736.svg)](https://zenodo.org/badge/latestdoi/373492736)

Script made in Matlab originally to extract the slope on FURA-2 experiments.

The data structure is a text file containing a different timepoint for each row and a different cell for each column.
The first column contains the time values and no x or y titles should be included in the text file.
There can be empty columns, they will be treated as empty data and the slope will be "NaN".
An important consideration regarding the number of row that need to be the same for each column.
You can duplicate the last value of your column until the last row to go around this limitation (using Excel for example).
Decimal numbers (floats) need to use the dot as separator and not the comma.

The script will open all the traces on a graph (median filter is applied and can be adjusted on line 94 "MedLen = 5") and you can select the time range in which you would like to find a slope.
You can tick the box in case you would like to fit a negative regression to your data.
The algorithm is based on the first derivative of the traces from which the fastest part is then isolated to fit a linear regression and extract the slope "a" from the equation "y = ax +b". 
There is a pause between each trace to let you observe that the fit is correct, you can modify it on line 55, "pause(1)".

Video example: [Slope From Linear Regression](https://youtu.be/a8gxQElWiVE)
