#!/usr/bin/env python3.4
import sys,re

# show dictionary
def ppdic(d):
    for k in sorted(d.keys()):
        print("%s : %s" % (k, d[k]))

def pp_sigma(g, sigma, all_u, pi):
    for i in sorted(sigma.keys()):
        pi_sum = sum(pi[h][0] for h in g.is2hs[i])
        if pi_sum == 0 : 
            print("%s ( - ) : %s" % (i, sigma[i]))
        else:
            allh_sum = float(sum(all_u[h] * pi[h][0] for h in g.is2hs[i]))
            print("%s ( %s ) : %s" % (i, (allh_sum/pi_sum), sigma[i]))

# cfr_plus
def makeSigma(q, isactions):
    sigma = {'' : {'-' : 1}}
    for i in isactions.keys():
        if  i == 'M,4 1 1 2,vR,gRR-,dR,vR,g.W+,dM,vM,gRR-,dC,vC':
            print("i=%s, s=%s" % (i,s))
        s = sum(max(0, v) for v in q[i].values())
        if  i == 'M,4 1 1 2,vR,gRR-,dR,vR,g.W+,dM,vM,gRR-,dC,vC':
            print("i=%s, s=%s" % (i,s))
        if s == 0: 
            sigma[i] = {k : (1.0 / len(q[i])) for k in isactions[i] }
        else:
            sigma[i] = {a : max(0, q[i][a])/float(s) for a in q[i].keys()}
    return sigma

def piMul(p1, p2):
    return tuple(p1[i]*p2[i] for i in range(3))
def makePiURec(g, h, p, pi, all_u, sigma):
    pi[h] = p
    if g.isTerminal(h):
        all_u[h] = g.utility(h)
        return all_u[h]
    ns = g.nextWithProb(h, sigma)
    all_u[h] = sum(makePiURec(g, h1, piMul(p, p1), pi, all_u, sigma) * p1[0] for h1, p1 in ns)
    return all_u[h]

def makePiU(g, sigma, pi, all_u):
    all_u[g.init] = makePiURec(g, g.init, (1, 1, 1), pi, all_u, sigma)
    return (pi, all_u)

def makeRDiff(g, h, i, a, sigma, pi, all_u):
    old_u = all_u[h]
    both_i = [g.makeIS(h, p) for p in range(2)]
    b_both = [{a : 1} if g.i2p(i) == p else sigma[both_i[p]] for p in range(2)]
    new_u = sum(all_u[n1] * p1[0] for n1,p1 in g.next(h, b_both))
    r = pi[h][2-g.i2p(i)] * (new_u - old_u)
    return (r if g.i2p(i) == 0 else -r)

def oneIterationCFR(g, count, r, sum_pi, sum_sigma, sigma, all_u, debugprint):
    for p in range(1,-1,-1):
        sigma = makeSigma(r, g.isactions)
        pi = {}
        makePiU(g, sigma, pi, all_u)
        if p==0: # count > 0:
            for i, hs in g.is2hs.items():
                pi_i = sum(pi[h][0] for h in g.is2hs[i])
                sum_pi[i] += pi_i
                for a in g.isactions[i]:
                    sum_sigma[i][a] += pi_i * sigma[i][a]
        if debugprint:
            for h in sorted(g.allh):
                print("%s -> pi = %s all_u = %s" % (h, pi[h], all_u[h]))
        rdiff = {i : {a : 0 for a in actions} for i, actions in g.isactions.items()}
        if debugprint:
            print("count=%s" % count)
            print("sigma=")
            ppdic(sigma)
            print("g.isactions=%s" % g.isactions)
            if p == 0: # count > 0:
                av_sigma = {i : { a : 0 if sum_pi[i] == 0 else sum_sigma[i][a] / float(sum_pi[i]) for a in actions } for i, actions in g.isactions.items() }
                print("av_sigma=")
                ppdic(av_sigma)
            print("pi,all_u=")
        for i, hs in g.is2hs.items():
            if g.i2p(i) != p: continue
            if debugprint:
                print("%s :" % i)
            for h in g.is2hs[i]:
                if debugprint:
                    print("  %s : pi=%s all_u=%s" % (h, pi[h], all_u[h]))
                for a in g.isactions[i]:
                    rd = makeRDiff(g,h, i, a, sigma, pi, all_u)
                    if debugprint:
                        print("    %s -> %s :" % (a, rd))
                    r[i][a] += rd
    return sigma

def cfr(g, iteration=100, debugprint=True):
    r = {i : { a : 0 for a in actions } for i, actions in g.isactions.items() }
    sum_pi = {i : 0 for i in g.isactions.keys() }
    sum_sigma = {i : { a : 0 for a in actions } for i, actions in g.isactions.items() }
    sigma = {}
    all_u = {}
    for count in range(iteration):
        sigma = oneIterationCFR(g, count, r, sum_pi, sum_sigma, sigma, all_u, debugprint)
    av_sigma = {i : { a : sigma[i][a] if sum_pi[i] == 0 else sum_sigma[i][a] / float(sum_pi[i]) for a in actions } for i, actions in g.isactions.items() }
    print("sigma=")
    ppdic(sigma)
    print("av_sigma=")
    ppdic(av_sigma)

def make_av_sigma(sum_sigma):
    r = {}
    for i, v in sum_sigma.items():
        psum = sum(v.values())
        r[i] = {a : 0 if psum == 0 else v[a] / float(psum) for a in v}
    return r
def oneIteration(g, count, q, sum_pi, sum_sigma, sigma, all_u, debugprint):
    for p in range(1,-1,-1):
        sigma = makeSigma(q, g.isactions)
        pi = {}
        makePiU(g, sigma,pi,all_u)
        if True: # p == 0: # count > 0:
            for i, hs in g.is2hs.items():
                if g.i2p(i) == p:
                    continue
                pi_i = sum(pi[h][2 - p] for h in g.is2hs[i])
                for a in g.isactions[i]:
                    sum_sigma[i][a] += pi_i * sigma[i][a] * (count + 1)
        if debugprint:
            print("count=%s" % count)
            print("sigma=")
            ppdic(sigma)
            print("pi=%s, all_u=%s" %(pi, all_u))
            print("g.isactions=%s, sum_pi=%s" % (g.isactions, sum_pi))
            print("pi,all_u=")
            for h in sorted(g.allh):
                print("%s -> pi = %s all_u = %s" % (h, pi[h], all_u[h]))
        rdiff = {i : {a : 0 for a in actions} for i, actions in g.isactions.items()}
        for i, hs in g.is2hs.items():
            if g.i2p(i) != p: continue
            if debugprint: print("%s :" % i)
            for h in g.is2hs[i]:
                if debugprint:
                    print("  %s : pi=%s all_u=%s" % (h, pi[h], all_u[h]))
                rdsum = 0
                for a in g.isactions[i]:
                    rd = makeRDiff(g, h, i, a, sigma, pi, all_u)
                    if debugprint:
                        print("    %s -> %s :" % (a, rd))
                    rdiff[i][a] += rd
                    rdsum += rd
                if debugprint:
                    if abs(rdsum) > 1e-8:
                        for a in g.isactions[i]:
                            print("rdiff %s %s %s %s" %(i, h, a, makeRDiff(g,h, i, a, sigma, pi, all_u)))
        for i, actions in g.isactions.items():
            q[i] = {a : max(0, q[i][a] + rdiff[i][a]) for a in actions}
            if debugprint:
                print("r[%s] = %s" %(i, q[i]))
    if debugprint:
        print("av_sigma=%s" % make_av_sigma(sum_sigma))
    return sigma

def cfrplus(g, iteration = 100, debugprint = True):
    q = {i : { a : 0 for a in actions } for i, actions in g.isactions.items() }
    sum_pi = {i : 0 for i in g.isactions.keys() }
    sum_sigma = {i : { a : 0 for a in actions } for i, actions in g.isactions.items() }
    sigma = {}
    all_u = {}
    for count in range(iteration):
        sigma = oneIteration(g, count, q, sum_pi, sum_sigma, sigma, all_u, debugprint)
    av_sigma = make_av_sigma(sum_sigma)
    pi = {}
    makePiU(g, sigma, pi, all_u)
    print("game value = %s, win percent = %s" % (all_u[g.init], (all_u[g.init] + 1) * 100.0 / 2.0))
    print("av_sigma=")
    pp_sigma(g, av_sigma, all_u, pi)
