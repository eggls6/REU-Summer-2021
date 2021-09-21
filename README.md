# REU-Summer-2021
This project assesses the magnitude of obstruction of Solar System object observations by current Starlink satellites in the LSST field of view.

## Table of contents
* [Background](#background)
* [How to Use](#how-to-use)
* [Technologies](#technologies)
* [Attributions](#attributions)
* [License](#license)

## Background
Vera C. Rubin Observatory will be a key facility for small body science in planetary astronomy over the next decade. It will carry out the Legacy Survey of Space and Time (LSST), observing the sky repeatedly over the course of ten years using a 6.5 m effective diameter telescope with a 9.6 square degree field of view, reaching approximately r = 24.5 mag per visit. The resulting dataset will provide extraordinary opportunities for both discovery and characterization of large numbers (10–100 times more than currently known) of small solar system bodies, furthering studies of planetary formation and evolution. [1]

Existing and planned satellite constellations threaten the success of LSST, Starlink being the largest contributor. At a minimum, a fraction of the area being imaged is lost to the sun-illuminated trails or significantly reduced in S/N (signal-to-noise ratio). [2]

This project uses only the currently orbiting 1600 Starlink satellites currently in orbit; with another 12,000 already approved for launch and 30,000 more pending approval, the magnitude of interference determined by this project would only increase dramatically. 

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
[1] Jones, R. L., Bannister, M. T., Bolin, B. T., Chandler, C. O., Chesley, S. R., Eggl, S., Greenstreet, S., Holt, T. R., Hsieh, H. H., Ivezic, Z., Juric, M., Kelley, M. S., Knight, M. M., Malhotra, R., Oldroyd, W. J., Sarid, G., Schwamb, M. E., Snodgrass, C., Solontoi, M., &amp; Trilling, D. E. (2021). The Scientific Impact of the Vera C. Rubin Observatory’s Legacy Survey of Space and Time (LSST) for Solar System Science. Bulletin of the AAS, 53(4). [https://doi.org/10.3847/25c2cfeb.d8909f28](https://doi.org/10.3847/25c2cfeb.d8909f28)

[2] Constance	Walker, & Jeffrey	Hall, & Lori	Allen, & Richard	Green, & Patrick	Seitzer, & Tony	Tyson, & Amanda	Bauer, & Kelsie	Krafton, & James	Lowenthal, & Joel	Parriott, & Phil	Puxley, & Tim	Abbott, & Gaspar	Bakos, & John	Barentine, & Cees	Bassa, & Blakeslee, John & Andrew	Bradshaw, & Cooke, Jeff & Daniel	Devost, & Yoachim, Peter. (2020). Impact of Satellite Constellations on Optical Astronomy and Recommendations Toward Mitigations. Bulletin of the AAS. 52. [https://doi.org/10.3847/25c2cfeb.346793b8](https://doi.org/10.3847/25c2cfeb.346793b8). 
* The functions ```icrf2radec``` and ```radec2icrf``` were provided by [Dr. Siegfried Eggl](https://github.com/eggls6) as part of this project. 
* The TLE files for the satellites are created by [CelesTrak](https://celestrak.com/NORAD/elements/).
* The functions for calculating bearing and cross track distance were adapted from [Chris Veness](http://www.movable-type.co.uk/scripts/latlong.html), available under the [MIT license](https://opensource.org/licenses/MIT). 

## License

