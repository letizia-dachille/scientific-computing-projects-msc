""" 
Creator: Letizia D'Achille 

Functions:

	kmer_search(S, k, freq, x)
		For each frequent k-mer computes an ordered list of its locations in S and the total number of its point-x mutations
	spet_location(S, k, p)
		Computes a SPET location in S for the parameters p, k
	get_my_pqs()
		Returns three pairs (p, q), where q is a SPET location in S for the parameters p, k = 40
"""

def kmer_search(S, k, freq, x):
    """
	For each frequent k-mer computes an ordered list of its locations in S and the total number of its point-x mutations

		Parameters:
			S (string): An n-letter string made up of the letters "A", "C", "G", "T"
			k (int): An integer 0 < k < n
            freq (int): A non-negative integer
            x (int): An integer 0 <= x < n
		
		Returns:
			L1 (list of lists of int): List containing for each frequent k-mer an ordered list of its locations in S
            L2 (list of int): List containing, at position i, the total number of point-x mutations of the k-mer corresponding to the i-th entry in L1
	"""
    d = {}
    n = len(S)
    for i in range(n-k+1):
        if (d.get(S[i:i+k]) == None):
            d[S[i:i+k]] = [i]
        else:
            d[S[i:i+k]] += [i]
    L1 = []
    L2 = []
    for key in d:
        if len(d[key]) >= freq:
            L1.append(d[key])
            count = 0
            for c in ["A", "C", "G", "T"]:
                if key[x] != c:
                    key1 = key[:x] + c + key[x + 1:]
                    if (d.get(key1) != None):
                        count += len(d[key1])
            L2.append(count)
    return L1, L2

""" 
# Test with Rabin-Karp
base = 4
mod = 100049

def extend_by_one(h, c):
    return (base * h + ord(c)) % mod

def remove_left(h, c, bstar):
    return (h - bstar * ord(c)) % mod

def hash_value(s):
    h = 0
    for c in s: h = extend_by_one(h, c)
    return h 
"""

def spet_location(S, k, p):
    """
	Computes a SPET location in S for the parameters p, k

		Parameters:
			S (string): An n-letter string made up of the letters "A", "C", "G", "T"
			k (int): An integer 0 < k < n
            p (int): An integer 0 <= p < n
		
		Returns:
			q (int): An integer with the SPET location in S for the parameters p, k
	"""
    d = {}
    n = len(S)
    q = None
    if (p - k > 0):
        """ 
        # Test with Rabin-Karp
        for i in range(p - k - 1, max(-1, p - 2*k - 1), -1):
            h = hash_value(S[i:i+k])
            if (d.get(h) == None):
                d[h] = [[S[i:i+k], i, 1]]
            else:
                flag = 0
                for j in range(len(d[h])):
                    if d[h][j][0] == S[i:i+k]:
                        d[h][j][2] += 1
                        flag = 1
                        break
                if flag == 0:
                    d[h].append([S[i:i+k], i, 1])
        bstar = pow(base, k-1, mod)    
        h = hash_value(S[:k])  
        for i in range(max(0, p - 2*k)):
            if (d.get(h) != None):
                for j in range(len(d[h])):
                    if d[h][j][0] == S[i:i+k]:
                       d[h][j][2] += 1
                       break
            if i < p - 2*k - 1:
                h = remove_left(h, S[i], bstar)
                h = extend_by_one(h, S[i+k]) 
        h = hash_value(S[p - k:p]) 
        for i in range(p - k, n - k + 1):
            if (d.get(h) != None):
                for j in range(len(d[h])):
                    if d[h][j][0] == S[i:i+k]:
                       d[h][j][2] += 1
                       break
            if i < n - k:
                h = remove_left(h, S[i], bstar)
                h = extend_by_one(h, S[i+k]) 
        qfreq = None
        for _, list in d.items():
            for el in list:
                key = el[0]
                count = 0
                for c in key:
                    if c == "C" or c == "G":
                        count += 1
                if (el[2] <= 5 and 0.35 * k <= count <= 0.65 * k):
                    if (q == None or el[2] < qfreq):
                        qfreq = el[2]
                        q = el[1]
        """
        for i in range(p - k - 1, max(-1, p - 2*k - 1), -1):
            if (d.get(S[i:i+k]) == None):
                d[S[i:i+k]] = [i, 1]
            else:
                d[S[i:i+k]][1] += 1
        for i in range(max(0, p - 2*k)):
            if (d.get(S[i:i+k]) != None):
                d[S[i:i+k]][1] += 1
        for i in range(p - k, n-k+1):
            if (d.get(S[i:i+k]) != None):
                d[S[i:i+k]][1] += 1
        qfreq = None
        for key in d:
            count = 0
            for c in key:
                if c == "C" or c == "G":
                    count += 1
            if (d[key][1] <= 5 and 0.35 * k <= count <= 0.65 * k):
                if (q == None or d[key][1] < qfreq):
                    qfreq = d[key][1]
                    q = d[key][0]
    return q

def get_my_pqs():
    """
	Returns three pairs (p, q), where q is a SPET location in S for the parameters p, k = 40

		Returns:
			L (list of pairs of int): A list of three pre-computed pairs (p, q)
	"""
    # q1 = None, q2 = number, p3 = MMDDYY = 022400
    return [(87852, None), (130849, 130808), (22400, 22359)]


if __name__ == "__main__":

    """
    # --------------- Clean files ---------------

    # Chloroplast of the tomato plant
    S = ""
    with open("sequence.fasta", "r") as f:
        for line in f:
            if line[0] != ">" and line[0] != ";":
                line = line.strip()
                S += ''.join([c for c in line.strip() if c in 'ACGT'])
    with open("sequence.txt", "w") as f: f.write(S)

    # Escherichia coli bacteria
    S = ""
    with open("sequence1.fasta", "r") as f:
        for line in f:
            if line[0] != ">" and line[0] != ";":
                line = line.strip()
                S += ''.join([c for c in line.strip() if c in 'ACGT'])
    with open("sequence1.txt", "w") as f: f.write(S)

    # Brewer’s yeast chromosome 7
    S = ""
    with open("sequence2.fasta", "r") as f:
        for line in f:
            if line[0] != ">" and line[0] != ";":
                line = line.strip()
                S += ''.join([c for c in line.strip() if c in 'ACGT'])
    with open("sequence2.txt", "w") as f: f.write(S)

    # Tomato chromosome 6
    S = ""
    with open("sequence3.fasta", "r") as f:
        for line in f:
            if line[0] != ">" and line[0] != ";":
                line = line.strip()
                S += ''.join([c for c in line.strip() if c in 'ACGT'])
    with open("sequence3.txt", "w") as f: f.write(S)

    # Human chromosome 21
    S = ""
    with open("sequence4.fasta", "r") as f:
        for line in f:
            if line[0] != ">" and line[0] != ";":
                line = line.strip()
                S += ''.join([c for c in line.strip() if c in 'ACGT'])
    with open("sequence4.txt", "w") as f: f.write(S)
    """

    """
    # --------------- Verify get_my_pqs() ---------------

    # Chloroplast of the tomato plant
    S = ""
    with open("sequence.txt", "r") as f: S = f.read()
    
    k = 40
    # q = None
    p1 = 87852
    # q = number
    p2 = 130849
    # p = MMDDYY = 022400
    p3 = 22400
    list = [(p1, spet_location(S, k, p1)), (p2, spet_location(S, k, p2)), (p3, spet_location(S, k, p3))]
    print(list == get_my_pqs())
    """

    """
    # --------------- Extra ---------------

    k = 10
    r = 10

    # Chloroplast of the tomato plant
    S = ""
    with open("sequence.txt", "r") as f: S = f.read()
    n = len(S)

    p = 7
    count = 0
    for _ in range(r):
        p = (83*p + 3643) % n
        q = spet_location(S, k, p)
        if (q == None):
            count += 1
    print(count)

    # Escherichia coli bacteria
    S = ""
    with open("sequence1.txt", "r") as f: S = f.read()
    n = len(S)

    p = 7
    count = 0
    for _ in range(r):
        p = (83*p + 3643) % n
        q = spet_location(S, k, p)
        if (q == None):
            count += 1
    print(count)

    # Brewer’s yeast chromosome 7
    S = ""
    with open("sequence2.txt", "r") as f: S = f.read()
    n = len(S)

    p = 7
    count = 0
    for _ in range(r):
        p = (83*p + 3643) % n
        q = spet_location(S, k, p)
        if (q == None):
            count += 1
    print(count)

    # Tomato chromosome 6
    S = ""
    with open("sequence3.txt", "r") as f: S = f.read()
    n = len(S)

    p = 7
    count = 0
    for _ in range(r):
        p = (83*p + 3643) % n
        q = spet_location(S, k, p)
        if (q == None):
            count += 1
    print(count)

    # Human chromosome 21
    S = ""
    with open("sequence4.txt", "r") as f: S = f.read()
    n = len(S)

    p = 7
    count = 0
    for _ in range(r):
        p = (83*p + 3643) % n
        q = spet_location(S, k, p)
        if (q == None):
            count += 1
    print(count)
    """
    