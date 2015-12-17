#!/usr/bin/env python2.7
import sys, re, fileinput, mafia
abs2ls = {}
def abstraction(i):
    if i[0] == 'G': 
        return 'G,'+str(mafia.isg2hist(i[2:]))+str(i.count('vM'))
    else: 
        return 'M,'+str(mafia.ism2hist(i[2:]))

def processline(l):
    if l[0] != 'G' and l[0] != 'M':
        return
    m = re.match(r'^([^\(]+) \( ([^\)]+) \) :', l)
#    m = re.match(r'^([^\(]+)(.*)' ,l)
    if m:
        i = m.group(1)
        v = m.group(2)
        if v == '-':
            return
        a = abstraction(i)
        if a not in abs2ls:
            abs2ls[a] = (float(v),l)
        else:
            if abs(float(v) - abs2ls[a][0]) > 1e-2:
                print("%s %s" % (l, abs2ls[a][1]))
    else:
        raise "format error %s" % l

for l in fileinput.input():
    processline(l)
