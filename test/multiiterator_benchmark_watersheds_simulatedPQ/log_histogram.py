#!/usr/bin/env python

src="log_1.log.bz2"

import bz2
import struct
from pylab import *
from numpy import *

f = bz2.BZ2File(src)
f.readline() # skip __BEGIN__ marker
s_shape = f.readline().split(" ")
shape = int(s_shape[1][:-1]), int(s_shape[2][:-1]), int(s_shape[3])
print "SHAPE=", shape
count=0
while f.readline():
    count+=1
count-=1 # discount end marker
print "COUNT=", count #4096766
f.seek(0)
f.readline() # skip __BEGIN__ marker
f.readline() # skip SHAPE

queuesize = 0
sizes = []
for i in xrange(count):
    l = f.readline().split(" ")
    op = int(l[0] == "__push:")
    if (op):
        queuesize += 1
    else:
        queuesize -= 1
    sizes.append(queuesize)
    prio = int(l[1][:-1])
    id = int(l[2])
f.close()

sizes = array(sizes)
h = hist(sizes)
