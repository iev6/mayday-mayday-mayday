import numpy as np

##Change input method from Javascript

coordx1, coordy1 = [0,0]
coordx2, coordy2 = [0,90]

start_city = [coordx1, coordy1]
dest_city = [coordx2, coordy2]


#calculation of average distance
ang1=np.pi/180*np.asarray(start_city, dtype= float)
ang2=np.pi/180*np.asarray(dest_city, dtype= float)

delx=np.cos(ang2[0])*np.cos(ang2[1])-np.cos(ang1[0])*np.cos(ang1[1])
dely=np.cos(ang2[0])*np.sin(ang2[1])-np.cos(ang1[0])*np.sin(ang1[1])
delz=np.sin(ang2[0])-np.sin(ang1[0])

c=np.sqrt(delx**2+dely**2+delz**2)

deldel=2.*np.arcsin(c/2)
dist=6350.*deldel #km

print dist
print 6350.*np.pi/2

speed=900. #km/hr
avgtime=dist/speed*60.   #minutes

# New York-Seattle  Flight ID of 1 to 20 characters
# 07/1995           Flight date (MM/YYYY)  
# KJFK              ICAO code of origin airport  
# KSEA              ICAO code of destination airport 
# 2                Number of en route altitudes
# 29               Minutes climbing to 1st en route altitude
# 35000    99      1st en route altitude:  feet   minutes
# 39000   147      2nd en route altitude:  feet   minutes
# 17               Minutes descending to destination airport

date_part = "04/2017"  ##
target = open("Input_exe.txt", 'w')
target.write("TestFlight\n")
target.write("04/2017\n")
target.write("<aircode1>\n")
target.write("<aircode2>\n")
target.write("1\n")
target.write("28\n")
target.write("39000 %f\n",avgtime)
target.write("17\n")
target.write("-----------------------\n")
