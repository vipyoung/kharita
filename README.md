# kharita
Robust and online map inference from crowded gps data

## Input
The input is a csv file in the following format:
vehicule_id,timestamp,lon,lat,speed,angle

Vehicule_id: important to identify trajectories.
Timestamp: important to sort gps points within trajectories
timestamp: in the format: yyyy-mm-dd hh:mm:ss+03
angle: in 0-360 interval. Angle to the north. 

## Running Kharita
Kharita can be invoked from command line as follows:

python khrita.py -p data -f data_2015-10-01 -r 25 -s 10 -a 40

-p: the folder containing the input data
-f: the input file name without its extension.
-r: the radius (cr) in meters used for the clustering
-s: the densification distance (sr) in meters
-a: the angle heading tolerance in degrees (0-360)

## Output
The code will produce a txt file containing the edges of the generated graph. 
