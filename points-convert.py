import sys
import csv

data = sys.argv[1]
natoms = int(sys.argv[2])

with open(data) as f:
    content = f.readlines()

content = [x.strip() for x in content]
#print content

xlist = content[0::3]
ylist = content[1::3]
zlist = content[2::3]
xyzlist = zip(range(natoms),xlist,ylist,zlist)

with open(data + '.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(xyzlist)
