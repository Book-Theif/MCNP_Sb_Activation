c Simulation of radioactive xenon via neutron capture
1 1 -6.691 1 -2 -3 $ Antimony 1 gram - good
2 2 -2.84 10 -11 -12 #1 $ Aluminum cylinder - good  
3 6 -0.001225 -20 21 -22 23 -24 25 #6 $ air space - good
4 3 -1.04 -30 31 -32 33 -34 35 #3 #1 #2 #6 #7 #8 $borated polyethylene
5 4 -0.95 -40 41 -42 43 -44 45 #1 #2 #3 #4 #6 #7 #8 $high density polyethylene
6 0 -13 $ Neutron Sphere
7 6 -0.001225 -16 -44 24 #1 #2 $ Air inside the beam
8 5 -1.44 16 -17 -44 24 #7 $ PVC Beam
99 0 40:-41:42:-43:44:-45 $ Outside world - good

c Antimony volume
1 1 pz -0.181171
2 1 pz 0.181171
3 1 c/z 26.5 0 0.36234
c Cylinders
10 1 pz -0.8
11 1 pz 0.8
12 1 c/z 26.5 0 1.7
c Neutron Sphere
13 1 s 0 0 0 .2
c PVC Pipe
16 1 c/x 0 0 1.905 
17 1 c/x 0 0 2.667
c Borated Polyethylene
20 1 pz 30.5
21 1 pz -30.5
22 1 py 15.35
23 1 py -15.35
24 1 px 23.5
25 1 px -23.5
c High Density Polyethylene
30 1 pz 56.5
31 1 pz -56.5
32 1 py 25.35
33 1 py -25.35
34 1 px 53.5
35 1 px -45.5
c outside
40 pz 74.3
41 pz -74.3
42 py 36.75
43 py -36.75
44 px 76.5
45 px -72.5

mode n
imp:n 1 7r 0
tr1 0 0 0 
m1  51121 0.5721 & 
    51123 0.4279 $Antimony
m11 51121 1
m12 51123 1
m2  13027 0.934 &
    29063 0.0203478 &
    29065 0.0456522 $Aluminum
m20 13027 1
m21 29063 1
m22 29065 1
m3  1001.70c -0.14286 &
    6012.70c -0.85714 $ Polyethylene
m4  1001.70c  -0.13571 &
    6012.70c -0.81428 &
    5010 -0.01 &
    5011 -0.04 $Borated polyethylene
m5 1001.70c -0.04848 &
   6012.70c -0.38432 &
   17035.70c -0.431072 &
   17037.70c -0.136128 $PVC
m6  7014.70c  -0.7722  &
    7015.70c  -0.0028  &
    8016.70c  -0.21053 & 
    18040.70c -0.0128  &
    1001.70c  -0.00137 &
    6000.70c  -0.0003  $ Air   
sdef pos=0 0 0 erg=2.45 cel=6
f4:n (1)
fm4 (-0.5721 11 102) (-0.4279 12 102)
f14:n (2)
fm14 (-0.934 20 102) (-0.0203478 21 102) (-0.0456522 22 102)
e4 1e-6 0.001 20 200
e14 1e-6 0.001 20 200
nps 1e7
