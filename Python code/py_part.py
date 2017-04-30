import numpy as np
import pyproj

##Change input method from Javascript
numberofsteps = 100

coordx1, coordy1 = [77.2, 28.617]
coordx2, coordy2 = [-99.117, 19.4]

start_city = [coordx1, coordy1]
dest_city = [coordx2, coordy2]

#Geodisc
g = pyproj.Geod(ellps='WGS84')
(az12, az21, dist) = g.inv(coordx1, coordy1, coordx2, coordy2)
# calculate line string along path with segments <= 1 km
lonlats = g.npts(coordx1, coordy1, coordx2, coordy2, 1 + int(dist / 100000))

# npts doesn't include start/end points, so prepend/append them
lonlats.insert(0, (coordx1, coordy1))
lonlats.append((coordx2, coordy2))

number_of_steps = len(lonlats)
#
# long_path = np.linspace(coordx1, coordx2, number_of_steps)  #Change it to geodisc
# lat_path = np.linspace(coordy1, coordy2, number_of_steps)

date_part = "K, 11, 2017/04/29, H0, D4, P0, C4, S0"  ##

target = open("Input_exe.txt", 'w')
target.write("START-----------------------\n")

for i in lonlats:
    if (i[0] < 0):
        arr = "S, "+str('{:f}'.format(abs(i[0])))+", "
    else:
        arr = "N, "+str('{:f}'.format(abs(i[0])))+", "

    if (i[1] < 0):
        arr = arr+"W, "+str('{:f}'.format(abs(i[1])))+", "
    else:
        arr = arr+"E, "+str('{:f}'.format(abs(i[1])))+", "

    arr = arr + date_part
    target.write(arr)
    target.write('\n')

target.write("STOP-----------------------\n")
