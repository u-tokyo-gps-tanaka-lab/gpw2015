#!/usr/bin/env python3.4
import sys,re,cfrplus,string

# white guard -> S (shiro)
(R,D,G,M,W,C,S) = (0,1,2,3,4,5,6)
HI = { 'R' : R, 'D' : D, 'G' : G, 'M' : M, 'W' : W, 'C' : C, 'S' : S}

# update hist with hs
def updateHist(hist, hs):
    for cmd in hs:
        if cmd[0]=='v':
            hist[HI[cmd[1]]] -= 1
        elif cmd[0] == 'g' and cmd[3] == '+':
            hist[HI[cmd[2]]] -= 1
        elif cmd == 'gRR-':
            hist[R] -= 1
            hist[C] += 1
        elif cmd == 'dR':
            hist[R] -= 1
            hist[W] += 1
        elif cmd == 'dC':
            hist[C] -= 1
            hist[W] += 1
        elif cmd == 'dG':
            hist[G] -= 1
            hist[S] += 1
    return hist

# make hist from history
def h2hist(h):
    hs = h.split(',')
    hist = [int(x) for x in hs[0].split(' ')] + [0] * 3
    return updateHist(hist, hs[1:])

# utility for history h
# if h is terminal and win by residence -> +1
# if h is terminal and win by mafia -> -1
# otherwise -> 0
def utility(h):
    hist = h2hist(h)
    if hist[M] == 0: 
        return 1
    elif hist[M] * 2 >= sum(hist): 
        return -1
    else: return 0

# h is a terminal history
def isTerminal(h):
    return utility(h) != 0

# get list of next histories of history h
def nextStates(h):
    hist = h2hist(h)
    hs=h.split(',')
    i3 = len(hs) % 3
    if i3 == 1: # Vote turn
        if len(hs) > 1 and hs[-1] == 'dM': return [h+',vM']
        return [h + ',v' + m for m in ['R', 'G', 'M', 'C'] if hist[HI[m]] > 0]
    elif i3 == 2: # Guard turn
        gs = ['-'] if hist[G] + hist[S] == 0 else [m for m in ['R', 'D', 'M', 'W', 'C'] if hist[HI[m]] > 0]
        r = []
        for m1 in ['R', 'D', 'G', 'W', 'C', 'S']:
            if hist[HI[m1]] > 0:
                for m in gs:
                    if m1 != m: 
                        r += [h+',g'+m+m1+'+']
                    else:
                        if hist[HI[m1]] > 1:
                            r += [h+',g'+m+m1+'+']
                        r += [h+',g'+m+m1+'-']
        return r
    else: # Detective turn
        if hist[D] == 0: return [h + ',d-']
        r = []
        for m in ['R', 'G', 'M', 'C']:
            if hist[HI[m]] > 0:
                if m != 'G' or hist[S] == 0:
                    r += [h+',d'+m]
        if hs[-1][3] == '+' and hs[-1][2] != 'W' and hs[-1][2] != 'S': r += [h+',d-']
        return r
def nextStatesGuardTurnWithProb(h, p_g, p_m):
    hist = h2hist(h)
    hs = h.split(',')
    # must be guard turn
    if h[-2] != 'v': return []
    p_g_r = {}
    for k, p in p_g.items():
        if k == 'R':
            s = hist[R] + hist[M]
            if hist[R] > 0:
                p_g_r['R'] = p * (hist[R]/float(s))
            if hist[M] > 0:
                p_g_r['M'] = p * (hist[M]/float(s))
        else:
            p_g_r[k] = p
    p_m_r = {}
    for k, p in p_m.items():
        if k == 'R':
            s = hist[R] + hist[G]
            if hist[R] > 0:
                p_m_r['R'] = p * (hist[R]/float(s))
            if hist[G] > 0:
                p_m_r['G'] = p * (hist[G]/float(s))
        elif k == 'W':
            s = hist[W] + hist[S]
            if hist[W] > 0:
                p_m_r['W'] = p * (hist[W]/float(s))
            if hist[S] > 0:
                p_m_r['S'] = p * (hist[S]/float(s))
        else:
            p_m_r[k] = p
    r = []
    for g_k, g_p in p_g_r.items():
        for m_k, m_p in p_m_r.items():
            mulp = g_p * m_p
            if g_k != m_k: 
                r += [(h+',g'+g_k+m_k+'+', (mulp, mulp, mulp))]
                continue
            mp = 1.0 / hist[HI[g_k]]
            if mp != 1.0:
                r += [(h+',g'+g_k+m_k+'+', (mulp * (1-mp), g_p * (1-mp), m_p * (1-mp)))]
            r += [(h+',g'+g_k+m_k+'-', (mulp * mp, g_p * mp, m_p * mp))]
    return r

def nextStatesWithProb(h, sigma):
    hist = h2hist(h)
    hs=h.split(',')
    i3 = len(hs) % 3
    if i3 == 1: # Vote turn
        if len(hs) > 1 and hs[-1] == 'dM':
            return [(h+',vM',(1,1,1))]
        r = []
        s = sum(hist[HI[m]] for m in ['R', 'G', 'M', 'C'])
        for m in ['R', 'G', 'M', 'C']:
            n_m = hist[HI[m]]
            if n_m > 0: 
                p = n_m/float(s)
                r += [(h+',v'+m, (p, p, p))]
        return r
    elif i3 == 2: # Guard turn
        g_i = makeISg(h)
        m_i = makeISm(h)
        p_g = {'-' : 1} if g_i not in sigma else sigma[g_i]
        p_m = sigma[m_i]
        return [(n, p) for n,p in nextStatesGuardTurnWithProb(h, p_g, p_m)]
    else: # detective turn
        if hist[HI['D']] == 0: return [(h+',d-',(1,1,1))]
        r = []
        s = sum(hist[HI[m]] for m in ['R', 'G', 'M', 'C'])
        if hs[-1][3] == '+' and hs[-1][2] in ['R', 'G', 'C']:
            s += 1
            r += [(h+',d-', (1.0/s, 1.0/s, 1.0/s))]
        for m in ['R', 'G', 'M', 'C']:
            n_m = hist[HI[m]]
            if n_m > 0:
                if m != 'G' or hist[S] == 0:
                    p=n_m/float(s)
                    r += [(h + ',d' + m, (p, p, p)) ]
        return r

def makeAllHistory(init):
    visited=set()
    prev = [init]
    while len(prev) > 0:
        next=[]
        for n in prev:
            visited.add(n)
            if isTerminal(n): continue
            for n1 in nextStates(n):
                next += [n1]
        prev = next
    return visited

def makeISm(h):
    if isTerminal(h) or h[-2] != 'v': return ''
    h = re.sub(r',g.(.)\+', r',g.\1+', h)
    return 'M,'+h.translate(string.maketrans('GS','RW'))

def ism2hist(i):
    hs=i.split(',')
    hist = [int(x) for x in hs[0].split(' ')]+[0]*3
    hist[R] += hist[G]
    hist[G] = 0
    hist[W] += hist[S]
    hist[S] = 0
    m_hist = updateHist(hist, hs[1:])
    return m_hist
def ism2action(i):
    m_hist = ism2hist(i)
    return [m for m in ['R', 'D', 'W', 'C'] if m_hist[HI[m]] > 0]

def isg2hist(i):
    hs=i.split(',')
    hist = [int(x) for x in hs[0].split(' ')]+[0]*3
    g_hist = updateHist(hist, hs[1:])
    g_hist[R] += g_hist[M]
    g_hist[M] = 0
    return g_hist

def isg2action(i):
    g_hist = isg2hist(i)
    return [m for m in ['R', 'D', 'W', 'C'] if g_hist[HI[m]] > 0]

def is2action(i):
    if i[0] == 'G': return isg2action(i[2:])
    else: return ism2action(i[2:])

# make the IS of guardsfrom history
def makeISg(h):
    hist = h2hist(h)
    if hist[G] + hist[S] == 0 or isTerminal(h) or h[-2] != 'v': return ''
    h = h.replace(',vM',',vR')
    h = h.replace(',gM',',gR')
    h = h.replace('dM,vR','dM,vM')
    return 'G,'+h

class Mafia:
    def __init__(self, nR, nD, nG, nM):
        self.init = "%d %d %d %d" % (nR, nD, nG, nM)
        self.allh = makeAllHistory(self.init)
        self.is2hs = {}
        for h in self.allh:
            if isTerminal(h): continue
            for i in [makeISg(h), makeISm(h)]:
                if i == '': 
                    continue
                elif i not in self.is2hs:
                    self.is2hs[i] = [h]
                else: self.is2hs[i].append(h)
        self.isactions = { i : is2action(i) for i in self.is2hs.keys()}

    def utility(self, h):
        return utility(h)
    def isTerminal(self, h):
        return isTerminal(h)
    def makeIS(self, h, player):
        return makeISg(h) if player == 0 else makeISm(h)
    def i2p(self, i):
        return (0 if i[0] =='G' else 1)
    def next(self, h, b_both):
        return nextStatesGuardTurnWithProb(h, b_both[0], b_both[1])
    def nextWithProb(self, h, sigma):
        return nextStatesWithProb(h, sigma)

def main():
    if len(sys.argv) < 5:
        raise Exception('Needs four arguments : nR nD nG nM')
    nR, nD, nG, nM = [int(x) for x in sys.argv[1:5]]
    nIterations = 100 if len(sys.argv) < 6 else int(sys.argv[5])
    g = Mafia(nR, nD, nG, nM)
    cfrplus.cfrplus(g, nIterations, False)

if __name__ == '__main__':
    main()

        
                
