#!/usr/bin/env python2.7
import sys, re, fileinput, mafia
abs2ls = {}
def abstraction(i):
    if i[0] == 'G': 
        return 'G,'+str(mafia.isg2hist(i[2:]))+str(i.count('vM'))
    else: 
        return 'M,'+str(mafia.ism2hist(i[2:]))+str(i.count('vM'))

def processline(l):
    if l[0] != 'G':
        return
    m = re.match(r'^([^\(]+) \( ([^\)]+) \) : .*\'R\': ([^,\}]+)[,\}]', l)
#    m = re.match(r'^([^\(]+)(.*)' ,l)
    if m:
        i = m.group(1)
        v = m.group(2)
        rval = m.group(3)
        if float(rval) < 0.99:
            return
        hist = mafia.isg2hist(i[2:])
        if hist[mafia.R] != 1:
            return
        if v == '-':
            return
#        if abs(float(v)) > 0.999:
#            return
        print(l)

for l in fileinput.input():
    processline(l)
