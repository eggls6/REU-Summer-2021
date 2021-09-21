# REU-Summer-2021
Determing the degree of obstruction of Solar System object observations by current Starlink satellites in the LSST field of view.

## Table of contents
* [Background](#background)
* [How to Use](#how-to-use)
* [Technologies](#technologies)
* [Attributions](#attributions)
* [License](#license)

## Background
This project assesses the magnitude of obstruction of Solar System object observations by current Starlink satellites in the LSST field of view.

## How to Use
Certain files are necessary to run the code. The block of code that imports this data and compiles it into a dataframe is listed below:
```
path = "/data/projects/lsst/baseline_fbs_v1p7p1/"
dir_list = os.listdir(path)

dflist=[]

for d in dir_list:
    if(d[0:2]=='S1'):
        dflist.append(pd.read_hdf('/data/projects/lsst/baseline_fbs_v1p7p1/'+d+'/visit-0000000.h5'))
        dflist.append(pd.read_hdf('/data/projects/lsst/baseline_fbs_v1p7p1/'+d+'/visit-0010000.h5'))
```
This can be replaced with data that the user wishes to work with, as long as the output dataframe ```dflist``` has columns for ``` ObjID ```, ``` FieldID ```, ``` FieldMJD ```, ``` AstRA(deg) ```, and ``` AstDec(deg) ``` with applicable values. Below lists what each variable represents:
* ``` ObjID ```: The ID number that represents each unique object
* ``` FieldID ```: The ID number that represents each unique field image
* ``` FieldMJD ```: The Modified Julian Date that each field image was taken
* ``` AstRA(deg) ```: The Right Ascension, in decimal degrees, of each object
* ``` AstDec(deg) ```: The Declination, in decimal degrees, of each object

The file ```Starlink tle.txt``` is a text file containing the TLE sets of currently orbiting Starlink satellites, last updated August 3, 2021 and sourced from [CelesTrak](https://celestrak.com/NORAD/elements/). Current data may be found at the same link.

The project checks each field against each satellite, testing if any objects are within 2 arcseconds position-wise and 30 seconds time-wise of the satellite trail. If yes, this is considered to be an obstructed object observation in the field. Only the objects that meet these conditions are plotted with the satellites.
	
## Technologies
Project is created with:
* Python version: 3.6.13

## Attributions
* The functions ```icrf2radec``` and ```radec2icrf``` were provided by [Dr. Siegfried Eggl](https://github.com/eggls6) as part of this project. 
* The TLE files for the satellites are created by [CelesTrak](https://celestrak.com/NORAD/elements/).
* The functions for calculating bearing and cross track distance were adapted from [Chris Veness](http://www.movable-type.co.uk/scripts/latlong.html), available under the [MIT license](https://opensource.org/licenses/MIT). 

## License

