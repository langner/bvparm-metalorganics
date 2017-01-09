import csv
import itertools
import random
import string
import sys


# The bvparm subdirectory contains code for parsing published literature values
# of bond valence model parameters, as gathered by I. David Brown, and provides
# them in the form of a nested Python dictionary.
sys.path.append('bvparm')
from bvparm import BondValenceParameters
bvparms = BondValenceParameters()

# References are denoted by letters, like in the original CIF file, and for this script
# we only want to consider using three possible sources:
#   a - Brown and Altermatt, (1985), Acta Cryst. B41, 244-247 (empirical)
#   b - Brese and O'Keeffe, (1991), Acta Cryst. B47, 192-197 (extrapolated)
#   j - Liu and Thorp (1993) Inorg. Chem. 32 4102-4105
refs_to_consider = ('a', 'b', 'j')

# These are the only anions we will consider here and refine BV parameters for.
anion_valences = {  'N' : -3, 'P' : -3,
                    'O' : -2, 'S' : -2, 'Se' : -2,
                    'F' : -1, 'Cl' : -1, 'Br' : -1, 'I' : -1 }
#anion_selection = ['N', "O", 'F', 'P', 'S', 'Cl', 'Se', 'Br', 'I']
#anion_selection = ['N', "O", 'F', 'S', 'Cl', 'Br']
#anion_selection = ['N', 'O']
anion_selection = ["O"]
nanions = len(anion_selection)

# This helper function takes an atom name and extracts the element name from it,
# for example it will turn C1 into just C and Cl25 into just Cl.
label2element = lambda lbl: "".join(itertools.takewhile(str.isalpha, lbl))

PREFIX = "data-split"


def get_position_in_list(L, S):
    """Get the index in S of the first element in the intersection of L and S."""
    for s in S:
        if s in L:
            return L.index(s)

def csv2distances(fname, new=False):
    """Read a CSV file and construct a dictionary that groups according to cation."""

    distances = {}
    with open(fname) as csvfile:

        csvreader = csv.reader(csvfile)
        header = map(string.strip, csvreader.next())
        for line in csvreader:

            if len(line) != len(header):
                print "WARNING: line does not have %i columns, skipping..." % len(header)
                print "Line is:", line
                continue

            irefcode = get_position_in_list(header, ['refcode', 'csd_accession_code'])
            refcode = line[irefcode].strip()
            ication_label = header.index('lab1') or header.index('atomname_cation')
            cation_label = line[ication_label].strip()
            lbl = "%s_%s" % (refcode, cation_label)
            if not lbl in distances:
                distances[lbl] = {}
            site = distances[lbl]

            ianion_lbl = header.index('lab2') or header.index('atomname_anion')
            anion_lbl = line[ianion_lbl].strip()
            idistance = header.index('dist1') or header.index('distance')
            distance = float(line[idistance].strip())
            if anion_lbl in site:
                args = (refcode, cation_lbl, anion_lbl)
                try:
                    assert site[anion_lbl] == distance
                    print "WARNING: two equal distances in %s for %s-%s" % args
                except AssertionError:
                    print "WARNING: two different distances in %s for %s-%s, creating new atom..." % args
                    anion_label = anion_lbl+"_double"
            site[anion_lbl] = distance

    return distances

def distances2csv(dists, params, fname, labels, picks, degen):

    with open(fname, "w") as csvfile:

        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['refcode', 'pick', 'cation_lbl', 'ox. state', 'BVS', 'anion_lbl', 'd', 'r0', 'bv'])

        tosave = []
        for id in range(degeneracy):
            for ip in [ip for ip,p in enumerate(picks) if p == id]:
                site_lbl = labels[ip]
                lbl = labels[ip]
                cation = dists[lbl]
                refcode, cation_lbl = site_lbl.split("_")
                av = 0.0
                bvs = []
                R = []
                for anion_lbl, d in cation.items():
                    cation_element = label2element(cation_lbl)
                    anion_element = label2element(anion_lbl)
                    r0 = params[anion_selection.index(anion_element) + id*len(anion_selection)]
                    R.append(r0)
                    bv = numpy.exp((r0-d)/0.37)
                    av += bv
                    bvs.append(bv)
                    #csvwriter.writerow([refcode.rjust(10), cation_lbl.rjust(10), " %i" % cation_valence, anion_lbl.rjust(10), " %.3f"%d, " %.3f" %bv, id, "%.3f"%r0])
                i = 0
                for anion_lbl, d in cation.items():
                    tosave.append([refcode, id, cation_lbl, cation_valence, av, anion_lbl, d, R[i], bvs[i]])
                    i += 1
        def bvcmp(a,b):
            if a[1] != b[1]:
                return cmp(int(a[1]), int(b[1]))
            elif a[4] != b[4]:
                return cmp(float(a[4]), float(b[4]))
            else:
                return cmp(float(a[6]), float(b[6]))
        tosave.sort(cmp=bvcmp)
        csvwriter.writerows(tosave)

def filter_sites(sites, cation_valence, bvparms=bvparms, anion_selection=anion_selection, valence_cutoff = 0.05):

    filtered = {}
    nsingle = 0
    nmissing = 0
    for site_lbl, site in sites.items():

        filtered[site_lbl] = {}

        refcode, cation_lbl = site_lbl.split("_")

        throw_away = False

        for anion_lbl, anion_dist in site.items():

            if throw_away:
                continue

            cation_element = label2element(cation_lbl)
            anion_element = label2element(anion_lbl)

            try:
                choices = bvparms[cation_element][cation_valence][anion_element][anion_valences[anion_element]]
                choose = [st for st in choices if st['ref'] in refs_to_consider]
            except KeyError:
                print "WARNING: in %s no appropriate BV parameters for %s-%s, skipping site..." % (site_lbl, cation_element, anion_element)
                throw_away = True
                continue

            r0 = choose[-1]['r0']
            b = 0.37
            d_cutoff = r0 - b*numpy.log(valence_cutoff*cation_valence)
            if anion_dist > d_cutoff:
                continue

            if anion_element not in anion_selection:
                print "WARNING: in %s, element %s found, throwing away..." % (site_lbl, anion_element)
                throw_away = True
                break

            filtered[site_lbl][anion_lbl] = anion_dist

        #if not "N" in [label2element(anion_lbl) for anion_lbl,dist in site.items()]:
        #    throw_away = True
        #    print 'WARNING: in %s no nitrogen, throwing away...' % site_lbl

        if len(filtered[site_lbl]) < 2:
            print "WARNING: only one anion for %s, throwing away..." % (site_lbl)
            nsingle += 1
            throw_away = True

        if throw_away:
            filtered.pop(site_lbl)

    print "Removed %i sites with only one anion" % nsingle
    return filtered


if __name__ == "__main__":

    import numpy
    from scipy import optimize

    import pylab
    plotting = "plot" in sys.argv

    cation_name = sys.argv[1]

    # This is normally 2 or 3, but 5 will denote we want to fit
    # both 2 and 3 together!
    cation_valence = int(sys.argv[2])

    # This label is only used in filenames below.
    cation_lbl = "%s%i" % (cation_name, cation_valence)

    degeneracy = 2
    nparams = nanions * degeneracy

    # The new files are created by Heping manually and have a slightly
    # different format than before, and contain only the column we need.
    #dists_all = csv2distances("initial_%s.csv" % cation_lbl)
    dists_all = csv2distances(PREFIX + "/bvparam2014d_%s.csv" % cation_lbl)
    
    dists = filter_sites(dists_all, cation_valence)
    picks = [random.choice(range(degeneracy)) for d in dists]

    print "Total number of sites:", len(dists_all)
    print "After processing:", len(dists)

    cation_element = label2element(dists.keys()[0].split("_")[1])
    start_bvparms = [bvparms[cation_element][cation_valence][a][anion_valences[a]] for a in anion_selection]
    start_bvparms = [[s for s in st if s['ref'] in refs_to_consider] for st in start_bvparms]
    start = [st[-1]['r0'] for id in range(degeneracy) for st in start_bvparms]

    el_in_site = [["".join(itertools.takewhile(str.isalpha, a)) for a in c.keys()] for c in dists.values()]
    element_counts = [sum([an in s for s in el_in_site]) for an in anion_selection]
    dist = dists

    anion_elements = anion_selection

    bond_valence = lambda r0, d, b: numpy.exp((r0 - d) / b)

    def site_valences(params, cation, ipick, elements, iparams):
        ipick_offset = ipick * nanions
        R0 = [params[ip + ipick_offset] for ip in iparams]
        b = 0.37
        #print R0, b
        bv = [bond_valence(R0[i], d, b) for i,d in enumerate(cation.values())]
        return bv

    labels = [lbl for lbl in dists.keys()]
    def valences(params, picks_local=picks):

        v = []
        r = []
        grad = [0.0 for p in params]
        bounds = [0.0 for p in params]

        for ilbl,lbl in enumerate(labels):
            cation = dists[lbl]
            elements = ["".join(itertools.takewhile(str.isalpha, anion_lbl)) for anion_lbl in cation.keys()]
            iparams = [anion_selection.index(el) for el in elements]
            ipick_offset = picks_local[ilbl] * nanions
            b = 0.37
            bv = site_valences(params, cation, picks_local[ilbl], elements, iparams)
            av = sum(bv)
            if not av:
                print av, bv, params
            for i,el in enumerate(elements):
                grad[iparams[i] + ipick_offset] += bv[i]*(av - cation_valence)
            for i,el in enumerate(set(elements)):
                bounds[iparams[i] + ipick_offset] += (b*numpy.log(cation_valence/av))**2

            v.append(av)
            r.append((av - 1.0*cation_valence)**2)

        b = 0.37
        grad = (2.0 / b) * numpy.array(grad)
        return numpy.array(v), numpy.array(r), grad, numpy.array(bounds)

    def grad(params):
        g = valences(params)[2]
#        if degeneracy > 1:
#            avg = 0.5 * (g[nanions+1:] + g[1:nanions])
#            g[nanions+1:] = avg
#            g[1:nanions] = avg
        return g

    def to_minimize(params):
        v, r, g, bounds = valences(params)
        return sum(r)

    import roman

    def make_initial_plot(fig, valences):
        fig.clf()
        ax = fig.add_subplot(111)
        h, b = numpy.histogram(valences, bins=100, range=[0.0, 2*cation_valence])
        width = b[1] - b[0]
        x = b[1:] - width/2.0
        ax.bar(x, h, width=width)
        ax.set_xlim([0.0, 2*cation_valence])
        ax.set_xlabel("bond valence sum $V_i$", fontsize=plot_labelsize)
        ax.set_ylabel("number of binding sites", fontsize=plot_labelsize)
        ax.text(1.5*cation_valence, 0.95*max(h), "%s(%s)" % (cation_element,roman.toRoman(cation_valence)), fontsize=32)
        return x, h

    def make_optimization_plot(fig, trends, valences):
        fig.clf()
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        for i in range(len(trends)):
            lbl = "%s-%s" % (cation_element, anion_selection[i % nanions])
            ax1.plot(range(len(trends[i])), trends[i], label=lbl)
        ax1.set_xlabel("optimization step", fontsize=plot_labelsize)
        ax1.set_ylabel("$R_0$ parameters", fontsize=plot_labelsize)
        ax1.legend(loc=6)
        h,b = numpy.histogram(valences, bins=100, range=[0.0, 2*cation_valence])
        width = b[1] - b[0]
        ax2.set_xlim([0.0, 2*cation_valence])
        ax2.set_xlabel("bond valence sum $V_i$", fontsize=plot_labelsize)
        ax2.set_ylabel("number of binding sites", fontsize=plot_labelsize)
        ax2.bar(b[1:] - width/2.0, h, width=width)
        return h,b

    plot_labelsize = 18
    if plotting:
        v, r, g, bounds = valences(start)
        fg1 = pylab.figure()
        x, h = make_initial_plot(fg1, v)
        s_anions = "-".join(anion_selection)
        fname = "plot_initial_%s_%s" % (cation_lbl, s_anions)
        if degeneracy > 1:
            fname += "_split"
        numpy.savetxt(fname + ".csv", zip(x, h), delimiter=',')
        fg1.savefig(fname + ".png")
        pylab.ion()
        pylab.show()

    # We do not want to do this anymore.
    #if cation_lbl == 'iron2':
    #    for lbl, cation in dist.items():
    #        av = 0.0
    #        for anion_lbl, d in cation.items():
    #            anion_element = "".join(itertools.takewhile(str.isalpha, anion_lbl))
    #            R0 = start[anion_selection.index(anion_element)]
    #            bv = numpy.exp((R0 - d)/0.37)
    #            av += bv
    #        if av > 3.0:
    #            dist.pop(lbl)
    #            labels.pop(labels.index(lbl))
    #            pass

    if plotting:
        to_plot = [[s] for s in start]
        fg2 = pylab.figure(figsize=(16,8))
        v, r, g, bounds = valences(start)
        x, h = make_optimization_plot(fg2, to_plot, v)
        pylab.draw()

    def ipick_dance(xk):
        flipped = 0
        for ilbl,lbl in enumerate(labels):
            cation = dists[lbl]
            elements = ["".join(itertools.takewhile(str.isalpha, anion_lbl)) for anion_lbl in cation.keys()]
            iparams = [anion_selection.index(el) for el in elements]
            ipick = picks[ilbl]
            av = sum(site_valences(xk, cation, ipick, elements, iparams))
            for id in range(1, degeneracy):
                ipick_alt = (ipick+id) % degeneracy
                av_alt = sum(site_valences(xk, cation, ipick_alt, elements, iparams))
                if abs(av-cation_valence) > abs(av_alt-cation_valence):
                    print "Switching pick (%i->%i) for %i (%s)" % (ipick, ipick_alt, ilbl, lbl)
                    picks[ilbl] = ipick_alt
                    flipped += 1
            if flipped > 0.1*len(labels):
                return

    def callback(xk):

        v, r, g, bounds = valences(xk)

        ipicks = [[i for i,ip in enumerate(picks) if ip == id] for id in range(degeneracy)]

        rpicks = [numpy.array([rr for irr,rr in enumerate(r) if irr in ip]) for ip in ipicks]
        ravg = [sum(rp)/len(rp) if len(rp)>0 else 0 for rp in rpicks]
        rargmax = [ipicks[id][rpicks[id].argmax()] if len(rpicks[id])>0 else -1 for id in range(degeneracy)]
        rmax =  [r[ip] if ip>=0 else 0.0 for ip in rargmax]

        gpicks = [numpy.array([gg for igg,gg in enumerate(g) if igg in ip]) for ip in ipicks]

        el_in_site = [[label2element(a) for a in c.keys()] for c in dists.values()]
        element_counts = numpy.array([sum([an in s for s in el_in_site]) for an in anion_selection])

        print ("%-7s " * nanions) % tuple(anion_selection)
        for id in range(degeneracy):
            x = xk[id*nanions:(id+1)*nanions]
            rmaxlbl = labels[rargmax[id]] if rargmax[id]>=0 else "-"
            print ("%-7.5f " * nanions) % tuple(x), "%4i" % len(ipicks[id]), "%.4f" % ravg[id], "%.2f (%s)" % (rmax[id], rmaxlbl), numpy.linalg.norm(gpicks[id])
        for id in range(degeneracy):
            bounds_ = bounds[id*nanions:(id+1)*nanions]
            b = numpy.sqrt(bounds_/element_counts)
            print ("%-7.5f " * nanions) % tuple(b)

        if plotting:
            for i in range(len(to_plot)):
                to_plot[i].append(xk[i])
            x, h = make_optimization_plot(fg2, to_plot, v)
            pylab.draw()

        ipick_dance(xk)

        cutoff = 3 #cation_valence
        while any(numpy.array(rmax) > cutoff**2):
            for id in range(degeneracy):
                if rmax[id] > cutoff**2:
                    ireject = rargmax[id]
                    to_reject = labels[ireject]
                    print "BVS deviation for %s is %.2f, more than %.2f, rejecting..." % (to_reject, rmax[id], cutoff)
                    dists.pop(to_reject)
                    picks.pop(ireject)
                    labels.pop(ireject)
                    r = numpy.delete(r, rargmax)
                    rargmax = [ipicks[id][rpicks[id].argmax()] if len(rpicks[id])>0 else -1 for id in range(degeneracy)]
                    rmax =  [r[ip] if ip>=0 else 0.0 for ip in rargmax]

    print "Homoleptic statistics"
    homoleptic = {}
    for lbl in dists:
        elements = ["".join(itertools.takewhile(str.isalpha, anion_lbl)) for anion_lbl in dists[lbl]]
        if len(set(elements)) == 1:
            el = elements[0]
            if not el in homoleptic:
                homoleptic[el] = []
            d = numpy.array(dists[lbl].values())
            b = 0.37
            homoleptic[el].append(b*numpy.log(1.0*cation_valence / sum(numpy.exp(-d/b))))
    for el,r in homoleptic.items():
        print el, len(r), numpy.average(r), numpy.std(r)

    print ("%-7s " * nanions) % tuple(anion_selection)

    print "***Intial callback***"
    callback(start)
    print "***Intial callback***"

    #opt = optimize.fmin_bfgs(to_minimize, start, fprime=grad, callback=callback, gtol=1e-06)
    opt = optimize.fmin(to_minimize, start, callback=callback, xtol=1e-06, ftol=1e-06)

    s_anions = "-".join(anion_selection)
    fname = 'optimized_%s_%s' % (cation_lbl, s_anions)
    if degeneracy > 1:
        fname += "_split"
    distances2csv(dists, opt, PREFIX + '/' + fname + ".csv", labels, picks, degeneracy)

    if plotting:
        fname = "plot_" + fname
        numpy.savetxt(PREFIX + '/' + fname + ".csv", zip(x, h), delimiter=',')
        fg2.savefig(PREFIX + '/' + fname + ".png")

    el_in_site = [["".join(itertools.takewhile(str.isalpha, a)) for a in c.keys()] for c in dists.values()]
    element_counts = []
    for id in range(degeneracy):
    	pick_el_in_site = [el_in_site[ip] for ip,p in enumerate(picks) if p==id]
    	element_counts.append([sum([an in s for s in pick_el_in_site]) for an in anion_selection])
    print "Final element counts:"
    for ec in element_counts:
        print zip(anion_selection, ec)

"""
homodists = {}

print "Homoleptic statistics"
homoleptic = {}
for lbl in dists:
    elements = ["".join(itertools.takewhile(str.isalpha, anion_lbl)) for anion_lbl in dists[lbl]]
    if len(set(elements)) == 1:
        homodists[lbl] = dists[lbl]
        el = elements[0]
        if not el in homoleptic:
            homoleptic[el] = []
        d = numpy.array(dists[lbl].values())
        b = 0.37
        homoleptic[el].append(b*numpy.log(cation_valence / sum(numpy.exp(-d/b))))
for el,r in homoleptic.items():
    print el, len(r), numpy.average(r), numpy.std(r)

print "********* HOMOLEPTIC OPTIMIZATION *********"

for lbl in dists.keys():
    if lbl not in homodists:
        dists.pop(lbl)
        labels.pop(labels.index(lbl))
opt = optimize.fmin_cg(to_minimize, start, fprime=grad, callback=callback, gtol=1e-03)
"""
