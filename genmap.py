'''
    Genetic Linkage Mapper
    Sergey Sysoev, 2013


'''
import sys
import itertools
import collections

def perr(*args):
    for x in args:
        print >>sys.stderr, x,
    print >>sys.stderr

OrganismRecord = collections.namedtuple("OrganismRecord",
    ["id", "parent1", "parent2", "sex", "allels"])

class Pedigree(object):
    class Organism(object):
        def __init__(self, record):
            self.id = record.id
            self.sex = record.sex
            self.allels = record.allels
            self.parents = []
            self.children = []

        def is_homozigota_at(self, i):
            return self.allels[2 * i] == self.allels[2 * i + 1]

        def __eq__(self, other):
            return isinstance(other, Pedigree.Organism) and self.id == other.id

        def __hash__(self):
            return hash(self.id)

    def __init__(self, M, number_of_species, locs_names, records):
        self.M = M
        self.number_of_species = number_of_species
        self.locs_names = locs_names
        organism_by_id = {r.id: Pedigree.Organism(r) for r in records}
        for r in records:
            child = organism_by_id[r.id]
            parents = [organism_by_id[p_id] for p_id in [r.parent1, r.parent2] if p_id]
            child.parents = parents
            for p in parents:
                p.children.append(child)
        self.organisms = organism_by_id.values()
        self.reveal_gametes()

    def reveal_gametes(self):
        def store_gamets(organism, g1, g2):
            organism.gamets1 = g1
            organism.gamets2 = g2

        # cycle through all species
        for s in self.organisms:
            # init gametes for the specie
            g1 = [0] * (self.M + 1)
            g2 = [0] * (self.M + 1)
            for i in range(self.M):
                if s.is_homozigota_at(i):
                    g1[i] = g2[i] = s.allels[2 * i]

            if not s.parents:
                store_gamets(s, g1, g2)
                continue

            # cycle through all loci
            p = s.parents[0]
            for i in range(self.M):
                if g1[i] != 0:                  # we already know this locus
                    continue
                if p.is_homozigota_at(i) and p.allels[2 * i] != 0: # parent is homozygota, but current specie is not
                    g1[i] = p.allels[2 * i]
                    g1[self.M] = p.id                # the parent for the gamete
                    g2[i] = 3 - p.allels[2 * i]     # 2->1, 1->2
                    if len(s.parents) == 2:
                        g2[self.M] = s.parents[1].id   # the parent for the gamete
                    continue
                else:
                    if len(s.parents) == 2:       # look at the other parent
                        m = s.parents[1]
                        if m.is_homozigota_at(i) and m.allels[2 * i] != 0:
                            g2[i] = m.allels[2 * i]
                            g2[self.M] = m.id
                            g1[i] = 3 - m.allels[2 * i]
                            g1[self.M] = p.id

            #
            #  31.05.13  Sysoev. It is useful to assign equal gametes to the parents, because homozygota child of
            #  heterozygota parent can be useful for further data retrival
            #
            if g1[self.M] == 0:
                g1[self.M] = s.parents[0].id
            if g2[self.M] == 0 and len(s.parents) == 2:
                g2[self.M] = s.parents[1].id

            store_gamets(s, g1, g2)

    def get_cistrans_matrix(self, stat=True, order_hint=None):
        # init the CIS - TRANS estimation matrix
        def zero_cistrans():
            return [[None] * self.M for _ in range(self.M)]

        def add_cistrans(a, b):
            def f(x):
                return (0, 0) if x is None else x
            return [[(ra + rb, nra + nrb)
                     for ((ra, nra), (rb, nrb)) in zip(map(f, ra), map(f, rb))]
                    for (ra, rb) in zip(a, b)]

        def organism_cistrans(o):
            het_loci = [i for i in range(self.M) if not o.is_homozigota_at(i)]
            ret = zero_cistrans()
            for (i, j) in itertools.product(het_loci, repeat=2):
                type1 = type2 = rec = nonrec = 0
                for ch in o.children:
                    gamete = {ch.gamets1[self.M]: ch.gamets1,
                              ch.gamets2[self.M]: ch.gamets2}[o.id]
                    if gamete[i] == 0 or gamete[j] == 0:    # the meiosis was uninformative on these loci
                        continue

                    if gamete[i] == gamete[j]:
                        type1 += 1                # AB or ab
                    else:
                        type2 += 1                # Ab or aB

                    # gather the reliable info
                    if (gamete[i], gamete[j]) in [(o.gamets1[i], o.gamets2[j]),
                                                  (o.gamets2[i], o.gamets1[j])]:
                        rec += 1                # RECOMBINATION
                    else:
                        nonrec += 1             # NO RECOMBINATION

                if stat and len(o.children) > 4:
                    if rec < min(type1, type2):
                        # print rec, nonrec, type1, type2
                        nonrec = max(type1, type2)
                        rec = min(type1, type2)

                ret[i][j] = (rec, nonrec)
            if order_hint:
                for i, left in enumerate(order_hint):
                    for j, right in enumerate(order_hint[i+1:], i+1):
                        prev = order_hint[j-1]
                        if ret[left][right] is None:
                            ret[left][right] = ret[right][left] = ret[left][prev]

            return ret

        cistranses = (organism_cistrans(o) for o in self.organisms if o.children)
        return reduce(add_cistrans, cistranses)

    #
    #   Given the recombinations, calculate the fractions
    #
    def get_pairwise_recombination_distance_matrix(self, order_hint=None):
        matrix = self.get_cistrans_matrix(order_hint=order_hint)
        fracs = [[1.0 * rec / (rec + nonrec) for (rec, nonrec) in row]
                 for row in matrix]
        return fracs

def not_empty_lines(f):
    return itertools.ifilter(lambda x: x,
        itertools.imap(lambda x: x.strip(), f))
#
#    Open and parse the .GEN file
#
def open_file(name):
    with open(name) as f:
        lines = not_empty_lines(f)

        next(lines)                 # number of FAMILIES - not used, assumed 1
        M = int(next(lines))       # number of loci
        locs_names = next(lines).split()  # names of loci

        next(lines)                 # read family number
        number_of_species = int(next(lines))

        records = []               # data array for the whole pedigree
        for _ in range(number_of_species):                  # READ FAMILY
            id, p1, p2, sex = map(int, next(lines).split())     # read species
            allels = [int(a) for a in next(lines).split()]
            records.append(OrganismRecord(id, p1, p2, sex, allels))

        return Pedigree(M, number_of_species, locs_names, records)

#
#    Given the recombination fractions matrix, try to form the order
#
def form_cluster(M, matrix):
    best_cluster = None
    best_cluster_len = 0

    #
    # try to form cluster from every node, than choose best
    #
    for locus in range(M):
        temp = [True] * M
        cluster = []
        temp[locus] = False                # mark as used
        cluster.append(locus)

        # find the first neighbor
        neighbor = None
        dist = 1
        for i in range(M):
            if temp[i] and matrix[locus][i] < dist:
                neighbor = i
                dist = matrix[locus][i]
        if not neighbor:
            if len(cluster) > best_cluster_len:
                best_cluster_len = len(cluster)
                best_cluster = cluster
            continue

        cluster.append(neighbor)
        temp[neighbor] = False

        # find the second neighbor
        neighbor2 = None
        dist = 1
        for i in range(M):
            if temp[i] and matrix[locus][i] < dist:
                neighbor2 = i
                dist = matrix[locus][i]
        if not neighbor2:
            if len(cluster) > best_cluster_len:
                best_cluster_len = len(cluster)
                best_cluster = cluster
            continue
        total_length = matrix[neighbor][neighbor2]

        # build the chain
        while True:
            locus = neighbor
            neighbor = None
            dist = 1
            for i in range(M):
                if temp[i] and matrix[locus][i] < dist and matrix[locus][i] <= total_length:
                    neighbor = i
                    dist = matrix[locus][i]
            if not neighbor:
                if len(cluster) > best_cluster_len:
                    best_cluster_len = len(cluster)
                    best_cluster = cluster
                break

            cluster.append(neighbor)
            temp[neighbor] = False
            #
            # fixing the situation with equal fractions
            #
            for i in range(M):
                if temp[i] and matrix[locus][i] <= dist and matrix[locus][i] <= total_length:
                    cluster.append(i)
                    temp[i] = False
            total_length += dist
    return best_cluster

#
#  We've formed the cluster, but oops.. some loci are not in it. Inserting them...
#
def insert_locus(cluster, locus, matrix):
    # find the closest neighbor
    neighbor1 = None
    dist1 = 1
    for i in range(len(cluster)):
        if matrix[cluster[i]][locus] < dist1:
            dist1 = matrix[cluster[i]][locus]
            neighbor1 = i

    if neighbor1 == 0:                        # near the beginning
        if matrix[locus][cluster[1]] > matrix[cluster[0]][cluster[1]]:
            cluster.insert(0, locus)        # our locus is the first
        else:
            cluster.insert(1, locus)        # no, it is the second
        return cluster
    if neighbor1 == len(cluster) - 1:        # near the end
        if matrix[locus][cluster[neighbor1 - 1]] > matrix[cluster[neighbor1]][cluster[neighbor1 - 1]]:
            cluster.insert(neighbor1 + 1, locus)        # it is last
        else:
            cluster.insert(neighbor1, locus)                # no, before the last
        return cluster

    # it is in the middle. But what is this neighbor? Right or Left?
    if matrix[locus][cluster[neighbor1 - 1]] < matrix[cluster[neighbor1]][cluster[neighbor1 - 1]]:
        cluster.insert(neighbor1, locus)
    else:
        cluster.insert(neighbor1 + 1, locus)
    return cluster

#
#    process the pedigree. Main function in the module
#         file_name - name of the CHR file with the pedigree data
#         order - order of loci that already known
#         stat - boolean value, whether to use statistical results
#
#    use it like:
#            process_pedigree("c:\\my_file.gen")
#       process_pedigree("c:\\my_file.gen", range(10), False) # first 10 loci are in the right order, use only the reliable results
#
def process_pedigree(file_name, order=None, stat=True):
    order = order or []
    pedigree = open_file(file_name)
    fracs = pedigree.get_pairwise_recombination_distance_matrix()
    for i in range(1):
        cluster = form_cluster(pedigree.M, fracs)
        for l in range(pedigree.M):
            if l not in cluster:
                cluster = insert_locus(cluster, l, fracs)
        fracs = pedigree.get_pairwise_recombination_distance_matrix(order_hint=cluster)


    for i in range(len(cluster)):
        name = pedigree.locs_names[cluster[i]]
        if i < len(cluster) - 1:
            print name, '   ', fracs[cluster[i]][cluster[i + 1]]
        else:
            print name

    return cluster, fracs


# NAME = 'c:\\python27\\Lib\\chr1000.gen'
NAME = 'chr1000.gen'
if len(sys.argv) > 1:
    NAME = sys.argv[1]


process_pedigree(NAME)
