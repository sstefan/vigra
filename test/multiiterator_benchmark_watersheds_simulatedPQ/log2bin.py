#!/usr/bin/env python

src="log_1.log.bz2"
dst="log.bin"

import bz2
import struct

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

f2 = file(dst, "wb")
f2.write(struct.pack("IIII",shape[0], shape[1], shape[2], count))
for i in xrange(count):
    l = f.readline().split(" ")
    op = int(l[0] == "__push:")
    prio = int(l[1][:-1])
    id = int(l[2])
    f2.write(struct.pack("III", op, prio, id))
f2.close()
f.close()
