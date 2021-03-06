#
#command line argument is name of a .inp file
#generate new .inp file with correct labels and units
#automatically generate a config file with default values
'''This file is part of Pipes.

    Pipes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Pipes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Pipes.  If not, see <http://www.gnu.org/licenses/>.
'''


import sys
import string
import math


def findit(keyword, lines):
    i1 = lines.index('['+keyword+']\n')
    where = [thing.find("\n") for thing in lines]
    i2 = where.index(0, i1+1)
    return (i1,i2)

# I change all c to 1 so that I do not need to convert units
def convert(x, units):
    c = 1.;
    if units =='f2m':     #feet to meters
        #c = 0.3048
	c=1
    elif units == 'm2f':  #meters to feet
        #c = 3.28084
	c=1
    elif units == 'fric': #convert friction factor
        #c = 0.01/140.
	c=1
    else:
        print 'no known conversion for' +units
    return [xi*c for xi in x]

def main(argv):
    f2m = 0.3048   #feet to meters
    dx = 1. #default grid spacing is 1 m
    a = 100#default pressure wave speed is 100 m/s
    T  =10. #detault run time is 1 s
    CFL = 0.8 #factor to make suggested value of M =number of time steps (hopefully) compliant with CFL
    
    s= sys.argv[1]
    print "Opening file %s" %s
    with open(s,'rU') as f:
        lines = f.readlines()
    (i1, i2) = findit('JUNCTIONS', lines) 
    orig = [it.split()[0] for it in lines[i1+2:i2]]        
    elevs= [float(it.split()[1]) for it in lines[i1+2:i2]]        
    (i1, i2) = findit('PUMPS', lines)
    orig += [it.split()[0] for it in lines[i1+2:i2]]
    elevs+= [float(it.split()[1]) for it in lines[i1+2:i2]]        
    (i1, i2) = findit('RESERVOIRS', lines)
    orig += [it.split()[0] for it in lines[i1+2:i2]]
    elevs+= [float(it.split()[1]) for it in lines[i1+2:i2]]        
    #need to keep track of reservoir head levels here and include them in bcs?
    print orig

    print elevs
    nodes = dict([(orig[i], str(i)) for i in range(len(orig))])
    (i3, i4)= findit('PIPES', lines)
    pipeinfo = [it.split()   for it in lines[i3+2:i4]]
    (i3, i4)= findit('COORDINATES', lines)
    xcoords = [float(it.split()[1]) for it in lines[i3+2:i4]]
    ycoords = [float(it.split()[2]) for it in lines[i3+2:i4]]

    #if coordinates are bullshit, scale them
    xmin = min(xcoords)
    xmax = max(xcoords)
    ymin = min(ycoords)
    ymax = max(ycoords)
    xscale = 100./(xmax-xmin)
    yscale = 100./(ymax-ymin)
    xcoords = [(x-xmin)*xscale for x in xcoords]
    ycoords = [(y-ymin)*yscale for y in ycoords]

    pipes = dict([(str(i),pipeinfo[i][0]) for i in range(len(pipeinfo))])
    conns = []
    ls = []
    ds = []
    rs = []

    for i in range(len(pipes)):
        conns.append(nodes[pipeinfo[i][1]])
        conns.append(nodes[pipeinfo[i][2]])
        ls.append(float(pipeinfo[i][3]))            #lengths
        ds.append(float(pipeinfo[i][4]))            #diameters
        rs.append(float(pipeinfo[i][5]))            #roughnesses
    print pipes
    print conns
    #assuming data is in GPM with length unit feet and roughness in H-W 
    #ds = [d/12 for d in ds] #assume D is in inches
    ds = [d for d in ds]	
    ls = convert(ls, 'f2m')
    ds = convert(ds, 'f2m')
    rs = convert(rs, 'fric')
    es = convert(elevs,'f2m')
    print "number of junctions is %d" %len(orig)
    print "number of pipes is %d" %len(pipeinfo)
    print len(elevs)
    print ds
    nodetypes = [conns.count(str(i)) for i in range(len(nodes))]
    print nodetypes
    q = s.find(".")
    newname = s[0:q]+'-2.0'+s[q:]
    fnew = open(newname, 'w')
    fnew.write('[TITLE]\n\n')
    fnew.write('[JUNCTIONS]\n;ID              	Elev        	Demand      	Pattern\n')
    for i in range(len(nodes)):
        fnew.write(' %d               	%3.2f         	0           	                	;\n'%(i, es[i]))
    fnew.write('\n[PIPES]\n;ID              	Node1           	Node2           	Length      	Diameter    	Roughness   	MinorLoss   	Status\n')    
    for i in range(len(pipes)):
        fnew.write(' %s               	  %s               	  %s               	%3.2f     	 %2.3f           %3.3f         	   0           	 Open  	;\n'  
                %(i, conns[2*i], conns[2*i+1], ls[i], ds[i], rs[i]))
    fnew.write('\n[COORDINATES]\n;Node            	X-Coord         	Y-Coord\n')
    for i in range(len(nodes)):
        fnew.write(' %s              	%5.2f              	%5.2f   \n'%(i, xcoords[i], ycoords[i]))
    
    fnew.close()
    fc = open(s[0:q]+'-2.0.config','w')
    fc.write('[PIPE_INFO]\n;-------------------------------------------\n;             initial  initial\n; ID    N       h       Q\n;-------------------------------------------\n')
    for i in range(len(pipes)):
        fc.write('%d      %d	  %2.3f       %1.f\n'%(i,math.ceil(ls[i]/dx), ds[i], 0))
    fc.write('[JUNCTION_INFO]\n;-------------------------------------------------------------------------------------\n;{-----for junction1s-----} | {--for junction2s--}             | {------for junction3s-------}|\n; ID   	jtype	bvaltype  bval   reflect   | offset   valveopen   | offset01   offset02   offset12 |\n;------------------------------------------------------------------------------------\n')
    for i in range(len(nodes)):       
        fc.write( '%2d        %s      1	     0       1	      0        1          0	     0		0    \n'%(i, nodetypes[i]))
    fc.write('\n[TIME_INFO];---------------------------------------\n;T (s)           M        Mi     a  (m/s)\n')
    M = int(a*T/(dx*CFL))
    fc.write(';-----------------------------------------\n%f         %d       10        %.2f'%(T,M,a))
    fc.close()
    print es
if __name__=="__main__":
    main(sys.argv[1:])
