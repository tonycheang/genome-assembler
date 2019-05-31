import sys
from array import array
from functools import lru_cache

''' Further testing of this data structure is required. '''
class CountMinSketch:
    # Various arrays of primes for testing purposes.
    # 20 primes cenetered around 6*10**7
    primes_6_10_7 = [59999879, 59999881, 59999887, 59999917, 59999957,
                     59999971, 59999981, 59999983, 59999993, 59999999,
                     60000011, 60000013, 60000023, 60000047, 60000049,
                     60000067, 60000071, 60000091, 60000103, 60000113]

    # 20 primes cenetered around 10**7
    primes_1_10_7 = [9999889, 9999901, 9999907, 9999929, 9999931,
                     9999937, 9999943, 9999971, 9999973, 9999991,
                     10000019, 10000079, 10000103, 10000121, 10000139,
                     10000141, 10000169, 10000189, 10000223, 10000229]

    # 20 primes cenetered around 5*10**6
    primes_5_10_6 = [4999871, 4999879, 4999889, 4999913, 4999933,
                     4999949, 4999957, 4999961, 4999963, 4999999,
                     5000011, 5000077, 5000081, 5000087, 5000101,
                     5000111, 5000113, 5000153, 5000161, 5000167]

    def __init__(self, num_rows):
        assert num_rows < len(
            CountMinSketch.primes_1_10_7), "Requested number of rows in CountMinSketch exceeds number of primes in list"
        self.num_rows = num_rows
        # Store in arrays with short ints to save space
        self.hash_values = [array('H', [0] * CountMinSketch.primes_1_10_7[i])
                            for i in range(num_rows)]

    def update(self, string, amount):
        string_hash = self._hash(string)
        for row in range(self.num_rows):
            self.hash_values[row][string_hash %
                                  len(self.hash_values[row])] += amount  # * self._get_sign(string_hash, row)

    def estimate(self, string):
        string_hash = self._hash(string)
        est_from_min = min(self.hash_values[row][string_hash %
                                                 len(self.hash_values[row])] for row in range(self.num_rows))
        return est_from_min

    @staticmethod
    @lru_cache(maxsize=512)
    def _hash(data, seed=0):
        ''' MurmurHash3 x86_32 port from Java by Maurus Decimus via StackOverflow'''
        c1 = 0xcc9e2d51
        c2 = 0x1b873593

        length = len(data)
        h1 = seed
        roundedEnd = (length & 0xfffffffc)  # round down to 4 byte block
        for i in range(0, roundedEnd, 4):
            # little endian load order
            k1 = (ord(data[i]) & 0xff) | ((ord(data[i + 1]) & 0xff) << 8) | \
                ((ord(data[i + 2]) & 0xff) << 16) | (ord(data[i + 3]) << 24)
            k1 *= c1
            k1 = (k1 << 15) | ((k1 & 0xffffffff) >> 17)  # ROTL32(k1,15)
            k1 *= c2

            h1 ^= k1
            h1 = (h1 << 13) | ((h1 & 0xffffffff) >> 19)  # ROTL32(h1,13)
            h1 = h1 * 5 + 0xe6546b64

        # tail
        k1 = 0

        val = length & 0x03
        if val == 3:
            k1 = (ord(data[roundedEnd + 2]) & 0xff) << 16
        # fallthrough
        if val in [2, 3]:
            k1 |= (ord(data[roundedEnd + 1]) & 0xff) << 8
        # fallthrough
        if val in [1, 2, 3]:
            k1 |= ord(data[roundedEnd]) & 0xff
            k1 *= c1
            k1 = (k1 << 15) | ((k1 & 0xffffffff) >> 17)  # ROTL32(k1,15)
            k1 *= c2
            h1 ^= k1

        # finalization
        h1 ^= length

        # fmix(h1)
        h1 ^= ((h1 & 0xffffffff) >> 16)
        h1 *= 0x85ebca6b
        h1 ^= ((h1 & 0xffffffff) >> 13)
        h1 *= 0xc2b2ae35
        h1 ^= ((h1 & 0xffffffff) >> 16)

        return h1 & 0xffffffff

    def __getitem__(self, key):
        # To allow working similar to dict
        return self.estimate(key)

    def __sizeof__(self):
        if hasattr(self, "total_mem"):
            return self.total_mem
        total_mem = 0
        total_mem += sys.getsizeof(self.num_rows)
        total_mem += sys.getsizeof(self.hash_values)
        for row in self.hash_values:
            total_mem += row.buffer_info()[1] * row.itemsize
        self.total_mem = total_mem
        self.total_mem += sys.getsizeof(self.total_mem)
        return total_mem
