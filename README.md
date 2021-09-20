# REU-Summer-2021
Determing the degree of obstruction of Solar System object observations by current Starlink satellites in the LSST field of view.

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Attributions](#attributions)

## General info
This project assesses the magnitude of obstruction of Solar System object observations by current Starlink satellites in the LSST field of view. Certain files are necessary to run the code. The block of code that imports this data and compiles it into a dataframe is listed below:
```
path = "/data/projects/lsst/baseline_fbs_v1p7p1/"
dir_list = os.listdir(path)

dflist=[]

for d in dir_list:
    if(d[0:2]=='S1'):
        dflist.append(pd.read_hdf('/data/projects/lsst/baseline_fbs_v1p7p1/'+d+'/visit-0000000.h5'))
        dflist.append(pd.read_hdf('/data/projects/lsst/baseline_fbs_v1p7p1/'+d+'/visit-0010000.h5'))
```
This can be replaced with data that the user wishes to work with, as long as the output dataframe ```dflist``` has columns for ``` ObjID ```, ``` FieldID ```, ``` FieldMJD ```, ``` AstRA(deg) ```, and ``` AstDec(deg) ``` with applicable values. Below lists what each variable represents.
* ``` ObjID ```: The ID number that represents each unique object
* ``` FieldID ```: The ID number that represents each unique field image
* ``` FieldMJD ```: The Modified Julian Date that each field image was taken
* ``` AstRA(deg) ```: The Right Ascension, in decimal degrees, of each object
* ``` AstDec(deg) ```: The Declination, in decimal degrees, of each object

The file ```Starlink tle.txt``` is a text file containing the TLE sets of currently orbiting Starlink satellites, last updated August 3, 2021 and sourced from [Celestrak](https://celestrak.com/NORAD/elements/). Current data may be found at the same link.
	
## Technologies
Project is created with:
* Python version: 3.6.13

## Attributions
* The functions ```icrf2radec``` and ```radec2icrf``` were provided by Dr. Siegfried Eggl as part of this project. 
* The functions for calculating bearing and cross track distance were 
