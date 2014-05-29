'''
    Genetic Linkage Mapper
    Sergey Sysoev, 2013


'''

NAME = 'chr1000.gen'

def read_line(f):
    for l in f:
        if l != '\n':
            return l.strip()

#
#    Open and parse the .GEN file
#
def open_file(name):
    try:
        f = open(name, 'r')
        read_line(f)                # number of FAMILIES - not used, assumed 1
        M = int(read_line(f))       # number of loci
        locs = read_line(f).split(' ')  # names of loci

        pedigree = []               # data array for the whole pedigree
        read_line(f)                # read family number
        number = int(read_line(f))  # number of species

        for i in range(number):                  # READ FAMILY
            specie = read_line(f).split(' ')     # read species
            allels = list(int(a) for a in read_line(f).split(' '))
            specie.append(allels)
            pedigree.append(specie)

        f.close()
        return M, number, locs, pedigree
    except:
        if f:
            f.close()
        return 0, 0, None, None

#
#    Reveal the gametes
#
def reveal_gametes(M, number, pedigree):
    # cycle through all species
    for s in pedigree:
        # init gametes for the specie
        g1 = range(M + 1)
        g2 = range(M + 1)
        for i in range(M):
            if s[4][2 * i] == s[4][2 * i + 1]:      # homozygota, gametes are equal
                g1[i] = g2[i] = s[4][2 * i]
            else:
                g1[i] = g2[i] = 0
            g1[M] = g2[M] = 0                       # we don't know the parent for the gamete yet

        # find parents
        parents = []
        for p in pedigree:
            if p[0] == s[1] or p[0] == s[2]:
                parents.append(p)
        if len(parents) == 0:
            s.append(g1)
            s.append(g2)
            continue

        # cycle through all loci
        p = parents[0]
        for i in range(M):
            if g1[i] != 0:                  # we already know this locus
                continue
            if p[4][2 * i] == p[4][2 * i + 1] and p[4][2 * i] != 0: # parent is homozygota, but current specie is not
                g1[i] = p[4][2 * i]
                g1[M] = p[0]                # the parent for the gamete
                g2[i] = 3 - p[4][2 * i]     # 2->1, 1->2
                if len(parents) == 2:
                    g2[M] = parents[1][0]   # the parent for the gamete
                continue
            else:
                if len(parents) == 2:       # look at the other parent
                    m = parents[1]
                    if m[4][2 * i] == m[4][2 * i + 1] and m[4][2 * i] != 0:
                        g2[i] = m[4][2 * i]
                        g2[M] = m[0]
                        g1[i] = 3 - m[4][2 * i]
                        g1[M] = p[0]

	#
	#  31.05.13  Sysoev. It is useful to assign equal gametes to the parents, because homozygota child of
	#  heterozygota parent can be useful for further data retrival
	#
	if g1[M] == 0:
	    g1[M] = parents[0][0]
	if g2[M] == 0 and len(parents) == 2:
            g2[M] = parents[1][0]

	# store gamets
        s.append(g1)
        s.append(g2)
    return



#
#    Process the gametes
#	 M - number of loci
#	 pedigree - species data with the phase revealed
#	 stat - boolean flag. If true then we try to estimate recombination fraction statistically among the children of one specie
#
def process_gametes(M, pedigree, stat = True):
    # init the CIS - TRANS estimation matrix
    cistrans = []
    for i in range(M):
        row = list([0, 0] for j in range(M))
        cistrans.append(row)

    # cycle through the species
    for s in pedigree:
        # find children
        children = []
        for ch in pedigree:
            if ch[1] == s[0] or ch[2] == s[0]:
                children.append(ch)

        if len(children) == 0:              # no children - no meioses
            continue

        # cycle through all locus combinations
        for i in range(M):
            a1 = s[4][i * 2]
            a2 = s[4][i * 2 + 1]
            if a1 == a2:                    # homozygota
                continue
            for j in range(i + 1, M):
                b1 = s[4][j * 2]
                b2 = s[4][j * 2 + 1]
                if b1 == b2:                # homozygota
                    continue

		type1 = type2 = rec = nonrec = 0

                # cycle through the children
                for ch in children:
                    gamete = None
                    if ch[5][M] == s[0]:    # find our gamete in the child
                        gamete = ch[5]
                    if ch[6][M] == s[0]:
                        gamete = ch[6]
                    if gamete == None:      # cant determine which of gametes is ours. 31.05.13 Sysoev: now this should not happen!
                        continue
                    if gamete[i] == 0 or gamete[j] == 0:    # the meiosis was uninformative on these loci
                        continue

		    # gather the statistical info
                    if gamete[i] == gamete[j]:
			type1 += 1		# AB or ab
		    else:
			type2 += 1		# Ab or aB

		    # gather the reliable info
                    if (gamete[i] == s[5][i] and gamete[j] == s[6][j]) or (gamete[i] == s[6][i] and gamete[j] == s[5][j]):
                        rec += 1	        # RECOMBINATION
                    else:
                        nonrec += 1             # NO RECOMBINATION

		if stat == True and len(children) > 4:
                    if rec < min(type1, type2):
			# print rec, nonrec, type1, type2
			nonrec = max(type1, type2)
			rec = min(type1, type2)

                cistrans[i][j][0] += rec
                cistrans[i][j][1] += nonrec

    return cistrans

#
#    Given the recombination fractions matrix, try to form the order
#
def form_cluster(M, matrix):
    best_cluster = None
    best_cluster_len = 0

    #
    # try to form cluster from every node, than choose best
    #
    for l in range(M):
        temp = range(M)
        cluster = []
        locus = l
        temp[l] = -1		# mark as used
        cluster.append(locus)

        # find the first neighbor
        neighbor = None
        dist = 1
        for i in range(M):
            if temp[i] != -1 and matrix[locus][i] < dist:
                neighbor = i
                dist = matrix[locus][i]
        if not neighbor:
            if len(cluster) > best_cluster_len:
                best_cluster_len = len(cluster)
                best_cluster = cluster
            continue

        cluster.append(neighbor)
        temp[neighbor] = -1

        # find the second neighbor
        neighbor2 = None
        dist = 1
        for i in range(M):
            if temp[i] != -1 and matrix[locus][i] < dist:
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
                if temp[i] != -1 and matrix[locus][i] < dist and matrix[locus][i] <= total_length:
                    neighbor = i
                    dist = matrix[locus][i]
            if not neighbor:
                if len(cluster) > best_cluster_len:
                    best_cluster_len = len(cluster)
                    best_cluster = cluster
                break

            cluster.append(neighbor)
            temp[neighbor] = -1
            #
            # fixing the situation with equal fractions
            #
            for i in range(M):
                if temp[i] != -1 and matrix[locus][i] <= dist and matrix[locus][i] <= total_length:
                    cluster.append(i)
                    temp[i] = -1
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

    if neighbor1 == 0:			# near the beginning
        if matrix[locus][cluster[1]] > matrix[cluster[0]][cluster[1]]:
            cluster.insert(0, locus)	# our locus is the first
        else:
            cluster.insert(1, locus)	# no, it is the second
        return cluster
    if neighbor1 == len(cluster) - 1:	# near the end
        if matrix[locus][cluster[neighbor1 - 1]] > matrix[cluster[neighbor1]][cluster[neighbor1 - 1]]:
            cluster.insert(neighbor1 + 1, locus)	# it is last
        else:
            cluster.insert(neighbor1, locus)		# no, before the last
        return cluster

    # it is in the middle. But what is this neighbor? Right or Left?
    if matrix[locus][cluster[neighbor1 - 1]] < matrix[cluster[neighbor1]][cluster[neighbor1 - 1]]:
        cluster.insert(neighbor1, locus)
    else:
        cluster.insert(neighbor1 + 1, locus)
    return cluster

#
#   Given the recombinations, calculate the fractions
#
def process_matrix(M, matrix, number):
    # first of all, calculate the fractions
    fracs = []
    for i in range(M):
        row = []
        for j in range(M):
            if j == i:
                recombinants = 0
                nonrecombinants = 1
            else:
                if j < i + 1:		# initial matrix is triangle (see how we form it). Fracs matrix should be full
                    recombinants = matrix[j][i][0]
                    nonrecombinants = matrix[j][i][1]
                else:
                    recombinants = matrix[i][j][0]
                    nonrecombinants = matrix[i][j][1]

            if recombinants != 0 or nonrecombinants != 0:
                row.append(float(recombinants) / (recombinants + nonrecombinants))
            else:
                row.append(2)   # no data for these loci.
                print 'this happened! no data for loci ', i, j
        fracs.append(row)

    return fracs



#
#    process the pedigree. Main function in the module
#	 file_name - name of the CHR file with the pedigree data
#	 order - order of loci that already known
#	 stat - boolean value, whether to use statistical results
#
#    use it like:
#    	process_pedigree("c:\\my_file.gen")
#       process_pedigree("c:\\my_file.gen", range(10), False) # first 10 loci are in the right order, use only the reliable results
#
def process_pedigree(file_name, order = [], stat = True):
    M, number, locs, pedigree = open_file(file_name)
    reveal_gametes(M, number, pedigree)
    rec = process_gametes(M, pedigree, stat)
    fracs = process_matrix(M, rec, number)

    # get the cluster
    if len(order) == 0:
        cluster = form_cluster(M, fracs)
    else:
        cluster = order

    for l in range(M):
        try:
            index = cluster.index(l)
        except:
            cluster = insert_locus(cluster, l, fracs)

    for i in range(len(cluster)):
        if i < len(cluster) - 1:
            print locs[cluster[i]], '   ', fracs[cluster[i]][cluster[i + 1]]
        else:
            print locs[cluster[i]]

    return cluster, fracs

process_pedigree(NAME, [], False)
