import csv
import itertools
import string
import sys


sys.path.append('bvparm')
from bvparm import BondValenceParameters
bvparms = BondValenceParameters()

anion_selection = ['N', "O", 'F', 'P', 'S', 'Cl', 'Se', 'Br', 'I']
anion_valences = { 'N' : -3, 'O' : -2, 'F' : -1, 'P' : -3, 'S' : -2, 'Cl' : -1, 'Se' : -2, 'Br' : -1, 'I' : -1 }

label2element = lambda lbl: "".join(itertools.takewhile(str.isalpha, lbl))

def csv2distances(fname):

    distances = {}
    with open(fname) as csvfile:

        csvreader = csv.reader(csvfile)
        header = map(string.strip, csvreader.next())
        for line in csvreader:

            if len(line) != len(header):
                print "WARNING: line does not have %i columns, skipping..." % len(header)
                continue

            refcode = line[header.index('csd_accession_code')].strip()
            cation_label = line[header.index('atomname_cation')].strip()
            lbl = "%s_%s" % (refcode, cation_label)
            if not lbl in distances:
                distances[lbl] = {}
            site = distances[lbl]

            anion_lbl = line[header.index('atomname_anion')].strip()
            distance = float(line[header.index('distance')].strip())
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

def filter_sites(sites, cation_valence, bvparms=bvparms, anion_selection=anion_selection, valence_cutoff = 0.05):

    filtered = {}
    for site_lbl, site in sites.items():

        filtered[site_lbl] = {}

        refcode, cation_lbl = site_lbl.split("_")

        throw_away = False

        for anion_lbl, anion_dist in site.items():

            cation_element = label2element(cation_lbl)
            anion_element = label2element(anion_lbl)

            if anion_element not in anion_selection:
                print "WARNING: in %s, element %s found, throwing site away..." % (site_lbl, anion_element)
                throw_away = True
                break

            choices = bvparms[cation_element][cation_valence][anion_element][anion_valences[anion_element]]
            choose = [st for st in choices if st['ref'] in ('a', 'b', 'j')]
            if not choose:
                print "WARNING: no appropriate BV parameters for %s-%s..." % (cation_element, anion_element)
                sys.exit(1)
            r0 = choose[-1]['r0']
            b = 0.37
            d_cutoff = r0 - b*numpy.log(valence_cutoff*cation_valence)
            if anion_dist > d_cutoff:
                continue

            filtered[site_lbl][anion_lbl] = anion_dist

        if len(filtered[site_lbl]) < 2:
            print "WARNING: only one anion for %s, throwing away..." % (site_lbl)
            throw_away = True

        if throw_away:
            filtered.pop(site_lbl)

    return filtered


if __name__ == "__main__":

    import numpy
    from scipy import optimize

    import pylab
    plotting = "plot" in sys.argv

    cation_name = sys.argv[1]
    cation_valence = int(sys.argv[2])
    cation_lbl = "%s%i" % (cation_name, cation_valence)

    dists_all = csv2distances("distances_%s.csv" % cation_lbl)
    dists = filter_sites(dists_all, cation_valence)

    print "Total number of sites:", len(dists_all)
    print "After processing:", len(dists)

    cation_element = label2element(dists.keys()[0].split("_")[1])
    start = [bvparms[cation_element][cation_valence][a][anion_valences[a]] for a in anion_selection]
    start = [[s for s in st if s['ref'] in ('a', 'b', 'j')] for st in start]
    start = [st[-1]['r0'] for st in start]
    print start
    #sys.exit(1)

    el_in_site = [["".join(itertools.takewhile(str.isalpha, a)) for a in c.keys()] for c in dists.values()]
    element_counts = [sum([an in s for s in el_in_site]) for an in anion_selection]
    dist = dists

    labels = [lbl for lbl in dists.keys()]
    def valences(params, cation_valence=1.0*cation_valence):
        v = []
        r = []
        grad = [0.0]*len(params)
        bounds = [0.0]*len(params)
        for lbl in labels:
            cation = dists[lbl]
            elements = ["".join(itertools.takewhile(str.isalpha, anion_lbl)) for anion_lbl in cation.keys()]
            iparams = [anion_selection.index(el) for el in elements]
            R0 = [params[ip] for ip in iparams]
            b = 0.37
            bv = [numpy.exp((R0[i] - d)/0.37) for i,d in enumerate(cation.values())]
            av = sum(bv)
            for i,el in enumerate(elements):
                grad[iparams[i]] += bv[i]*(sum(bv) - cation_valence)
            for i,el in enumerate(set(elements)):
                bounds[iparams[i]] += (b*numpy.log(cation_valence/av))**2
            v.append(av)
            r.append((av - 1.0*cation_valence)**2)
        b = 0.37
        grad = (2.0 / b) * numpy.array(grad)
        return numpy.array(v), numpy.array(r), grad, numpy.array(bounds)

    def grad(params):
        return valences(params)[2]

    def to_minimize(params):
        v, r, g, bounds = valences(params)
        return sum(r)

    to_plot = [[s] for s in start]
    fg1 = pylab.figure(figsize=(16,8))
    ax1 = fg1.add_subplot(121)
    ax2 = fg1.add_subplot(122)
    for i in range(len(to_plot)):
        ax1.plot(range(len(to_plot[i])), to_plot[i])
    v, r, g, bounds = valences(start)
    h,b = numpy.histogram(v, bins=100, range=[0.0, 2*cation_valence])
    width = 0.8 * (b[1] - b[0])
    ax2.set_xlim([0.0, 2*cation_valence])
    ax2.set_xlabel("bond valence sum $V_i$")
    ax2.set_ylabel("count")
    ax2.bar(b[1:] + width/2.0, h, width=width)
    if plotting:
        pylab.ion()
        pylab.show()

    if cation_lbl == 'iron2':
        for lbl, cation in dist.items():
            av = 0.0
            for anion_lbl, d in cation.items():
                anion_element = "".join(itertools.takewhile(str.isalpha, anion_lbl))
                R0 = start[anion_selection.index(anion_element)]
                bv = numpy.exp((R0 - d)/0.37)
                av += bv
            if av > 3.0:
                dist.pop(lbl)
                labels.pop(labels.index(lbl))
                pass

    def callback(xk):
        v, r, g, bounds = valences(xk)
        ravg = sum(r)/len(r)
        rargmax = r.argmax()
        rmax =  r[rargmax]
        print "%-7.5f "*len(start) % tuple(xk), "%.4f" % ravg, "%.2f (%s)" % (rmax, labels[r.argmax()]), numpy.linalg.norm(g)
        print "%-7.5f "*len(start) % tuple(numpy.sqrt(bounds / element_counts))
        for i in range(len(to_plot)):
            to_plot[i].append(xk[i])
        fg1.clf()
        ax1 = fg1.add_subplot(121)
        ax2 = fg1.add_subplot(122)
        for i in range(len(to_plot)):
            ax1.plot(range(len(to_plot[i])), to_plot[i])
        h,b = numpy.histogram(v, bins=100, range=[0.0, 2*cation_valence])
        width = 0.8 * (b[1] - b[0])
        ax2.set_xlim([0.0, 2*cation_valence])
        ax2.set_xlabel("bond valence sum $V_i$")
        ax2.set_ylabel("count")
        ax2.bar(b[1:] + width/2.0, h, width=width)
        if plotting:
            pylab.draw()
        cutoff = cation_valence
        while rmax > cutoff**2:
            to_reject = labels[r.argmax()]
            print "BVS deviation for %s is %.2f, more than %.2f, rejecting..." % (to_reject, rmax, cutoff)
            dists.pop(to_reject)
            labels.pop(labels.index(to_reject))
            r = numpy.delete(r, rargmax)
            rargmax = r.argmax()
            rmax = r[rargmax]

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
            homoleptic[el].append(b*numpy.log(cation_valence / sum(numpy.exp(-d/b))))
    for el,r in homoleptic.items():
        print el, len(r), numpy.average(r), numpy.std(r)

    #print anion_elements
    print ("%-7s "*len(start)) % tuple(anion_selection)

    print "***Intial callback***"
    callback(start)
    print "***Intial callback***"
    opt = optimize.fmin_cg(to_minimize, start, fprime=grad, callback=callback, gtol=1e-03)

pylab.savefig("%s.optimized.png" % cation_name)

el_in_site = [["".join(itertools.takewhile(str.isalpha, a)) for a in c.keys()] for c in dists.values()]
element_counts = [sum([an in s for s in el_in_site]) for an in anion_elements]
print "Final element counts:", zip(anion_elements, element_counts)

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

for lbl in dists.keys():
    if lbl not in homodists:
        dists.pop(lbl)
        labels.pop(labels.index(lbl))
opt = optimize.fmin_cg(to_minimize, start, fprime=grad, callback=callback, gtol=1e-03)
"""