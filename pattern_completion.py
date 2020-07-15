import copy
import itertools as it
from collections import deque
import cProfile

class MPoly:
    """Multivariable polynomial"""
    def __init__(self, terms):
        if type(terms) is MPoly:
            self.terms = dict(terms.terms)
            self.len = terms.len
        else:
            if type(terms) is int:
                if terms == 0:
                    terms = {}
                else:
                    terms = {():terms}
            self.terms = terms
            for x in terms:
                self.len = len(x)
                break
            else:
                self.len = 0
    def __copy__(self):
        return MPoly(dict(self.terms))
    def coerce(self, other):
        if type(other) is not MPoly:
            other = MPoly(other)
        if self.len > other.len:
            return (self, other.augment_copy(self.len - other.len))
        if self.len < other.len:
            return (self.augment_copy(other.len - self.len), other)
        return (self, other)
    def __iter__(self):
        return iter(self.terms)
    def __add__(self, other):
        (self, other) = MPoly.coerce(self, other)
        out = dict(self.terms)
        for x,y in other.terms.items():
            if x in out:
                out[x] += y
                if out[x] == 0:
                    del out[x]
            else:
                out[x] = y
        return MPoly(out)
    def __radd__(self, other):
        return self + other
    def __sub__(self, other):
        (self, other) = MPoly.coerce(self, other)
        out = dict(self.terms)
        for x,y in other.terms.items():
            if x in out:
                out[x] -= y
                if out[x] == 0:
                    del out[x]
            else:
                out[x] = -y
        return MPoly(out)
    def __rsub__(self, other):
        return self - other
    # use this for adding lots of MPoly objects
    # to avoid unnecessary calls to dict.
    def sum(self, *other_iter):
        M = max(other_iter, key=lambda x:x.len, default=self)
        (self, _) = MPoly.coerce(self, M)
        out = dict(self.terms)
        for other in other_iter:
            (_, other) = MPoly.coerce(self, other)
            for x,y in other.terms.items():
                if x in out:
                    out[x] += y
                    if out[x] == 0:
                        del out[x]
                else:
                    out[x] = y
        return MPoly(out)
    def __mul__(self, other):
        (self, other) = MPoly.coerce(self, other)
        to_add = []
        for x,y in other.terms.items():
            to_add.append(MPoly({tuple(a+b for a,b in zip(z,x)):y*w
                                 for z,w in self.terms.items()}))
        if not to_add:
            return MPoly({})
        return MPoly.sum(*to_add)
    def __rmul__(self, other):
        return self * other
    def __neg__(self):
        return MPoly({x:-y for x,y in self.terms.items()})
    def __pow__(self, n): #using recursive repeated squaring
        if n == 0:
            return MPoly(1)
        if n == 1: #slight optimization
            return self
        if n < 0:
            raise ValueError("Negative not allowed for polynomial power.")
        if n != int(n):
            raise ValueError("Power must be integer.")
        if n % 2: #odd
            return self * (self * self) ** (n // 2)
        return (self * self) ** (n // 2) # even
    def __str__(self):
        if not self.terms:
            return "0"
        to_join = [("" if y == 1 else "-" if y == -1 else str(y))
                    + " ".join(f"a{i}^{j}" if j > 1
                               else f"a{i}"
                               for i,j in enumerate(x) if j >= 1)
                               for x,y in self.terms.items()]
        return " + ".join("1" if not x else "-1" if x == "-" else x for x in to_join)
    def __repr__(self):
        return str(self)
    def augment(self, n=1):
        self.terms = {x+(0,)*n:y for x,y in self.terms.items()}
        self.len += 1
    def augment_copy(self, n=1):
        return MPoly({x+(0,)*n:y for x,y in self.terms.items()})
    def solve(self):
        """Solve p > 0 and p < 0 in as many variables as possible."""
        if self.terms and all(x > 0 for x in self.terms.values()):
            return [True, False]
        if not self.terms or not any(x > 0 for x in self.terms.values()):
            return [False, True]
        not2 = [True]*self.len
        yes1 = [False]*self.len
        for t in self.terms:
            for i,x in enumerate(t):
                if x == 1:
                    yes1[i] = True
                elif x > 1:
                    not2[i] = False
        possible = [i for i,(a,b) in enumerate(zip(not2, yes1)) if a and b]
        out = [[], []]
        # a + b*x > 0 and a + b*x < 0
        for i in possible:
            a = MPoly({t:v for t,v in self.terms.items() if t[i] == 0})
            b = MPoly({tuple(0 if j == i else x for j,x in enumerate(t)):v
                       for t,v in self.terms.items() if t[i] == 1})
            apos = a.solve()
            bpos = b.solve()
            ## ... > 0
            # True if b > 0 and a > 0 (already dealt with)
            # False if b < 0 and a < 0 (already dealt with)
            # x > -a/b if b > 0 and a < 0
            # x < -a/b if b < 0 and a > 0
            if apos[0] and bpos[1]:
                out[0].append((i, '<', MRatio(a,-b), apos[0], bpos[1]))
            if apos[1] and bpos[0]:
                out[0].append((i, '>', MRatio(-a,b), apos[1], bpos[0]))
            ## ... < 0
            # False if b > 0 and a > 0 (already dealt with)
            # True if b < 0 and a < 0 (already dealt with)
            # x < -a/b if b > 0 and a < 0
            # x > -a/b if b < 0 and a > 0
            if apos[0] and bpos[1]:
                out[1].append((i, '>', MRatio(a,-b), apos[0], bpos[1]))
            if apos[1] and bpos[0]:
                out[1].append((i, '<', MRatio(-a,b), apos[1], bpos[0]))
        return tuple(out[0]), tuple(out[1])

class MRatio:
    """Multivariable rational function"""
    def __init__(self, num, den=1):
        if type(num) is MRatio:
            self.num = copy.copy(num.num)
            self.den = copy.copy(num.den)
        else:
            self.num = MPoly(num)
            self.den = MPoly(den)
        num_mins = []
        den_mins = []
        for i in range(min(self.num.len, self.den.len)):
            num_mins.append(min(x[i] for x in self.num.terms))
            den_mins.append(min(x[i] for x in self.den.terms))
        mins = [min(x,y) for x,y in zip(num_mins, den_mins)] + [0]*abs(self.num.len - self.den.len)
        self.num = MPoly({tuple(t[i] - mins[i] for i in range(len(t))):v for t,v in self.num.terms.items()})
        self.den = MPoly({tuple(t[i] - mins[i] for i in range(len(t))):v for t,v in self.den.terms.items()})
    def __add__(self, other):
        other = MRatio(other)
        return MRatio(self.num*other.den + self.den*other.num, self.den*other.den)
    def __radd__(self, other):
        return self + other
    def __sub__(self, other):
        other = MRatio(other)
        return MRatio(self.num*other.den - self.den*other.num, self.den*other.den)
    def __rsub__(self, other):
        return self - other
    def sum(self, *other_iter):
        return sum(other_iter, self)
    def __mul__(self, other):
        other = MRatio(other)
        return MRatio(self.num*other.num, self.den*other.den)
    def __rmul__(self, other):
        return self * other
    def __truediv__(self, other):
        other = MRatio(other)
        return MRatio(self.num*other.den, self.den*other.num)
    def __rtruediv__(self, other):
        return self / other
    def __neg__(self):
        return MRatio(-self.num, self.den)
    def __pow__(self, n): #using recursive repeated squaring
        if n == 0:
            return MRatio(1, 1)
        if n < 0:
            return MRatio(self.den ** n, self.num ** n)
        return MRatio(self.num ** n, self.den ** n)
    def __str__(self):
        return f"({self.num}) / ({self.den})"
    def __repr__(self):
        return str(self)
    def augment(self, n=1):
        self.num.augment(n)
        self.den.augment(n)
    def substitute(self, index, new):
        numden = []
        fast = 0
        for x in [self.num, self.den]:
            if index >= x.len:
                numden.append(MRatio(x))
            elif max(term[index] for term in x.terms) <= 1: #a + b * x
                a = MPoly({t:v for t,v in x.terms.items() if t[index] == 0})
                b = MPoly({tuple(0 if j == index else y for j,y in enumerate(t)):v
                           for t,v in x.terms.items() if t[index] == 1})
                numden.append(MRatio(a * new.den + b * new.num, new.den))
                fast += 1
            else:
                to_sum = []
                for t, c in x.terms.items():
                    pow_ = t[index] if index < len(t) else 0
                    quot = MRatio(MPoly({tuple(0 if i == index else y
                                               for i,y in enumerate(t)):c}), 1)
                    to_sum.append(quot * new ** pow_)
                numden.append(MRatio.sum(*to_sum))
        if fast == 2:
            return MRatio(numden[0].num, numden[1].num)
        return numden[0] / numden[1]
    def flipsign(self):
        self.num = -self.num
        self.den = -self.den

class Matrix:
    def __init__(self, data, depth=None):
        self.data = [[x if x is None else MRatio(x) for x in y] for y in data]
        self.total = max(x.num.len for y in self.data for x in y if x is not None)
        if depth is None:
            self.depth = self.total
        else:
            self.depth = depth
        self.len = (len(data), len(data[0]))
    def __len__(self):
        return len(self.data)
    def __iter__(self):
        return iter(self.data)
    def __getitem__(self, index):
        if type(index) is slice:
            return Matrix(self.data[index], self.depth)
        if type(index) is list:
            return Matrix([self.data[i] for i in index], self.depth)
        if type(index) is int:
            return self.data[index]
        if type(index[0]) is slice:
            if type(index[1]) is slice:
                return Matrix([x[index[1]] for x in self.data[index[0]]], self.depth)
            if type(index[1]) is list:
                return Matrix([[x[i] for i in index[1]] for x in self.data[index[0]]], self.depth)
            return [x[index[1]] for x in self.data[index[0]]]
        if type(index[0]) is list:
            if type(index[1]) is slice:
                return Matrix([self.data[i][index[1]] for i in index[0]], self.depth)
            if type(index[1]) is list:
                return Matrix([[self.data[i][j] for j in index[1]] for i in index[0]], self.depth)
            return [x[index[1]] for i in index[0] for x in self.data[i]]
        return self.data[index[0]][index[1]]
    def __str__(self):
        return '\n[' + ']\n ['.join('\t'.join(str(y) for y in x) for x in self.data) + ']]\n'
    def __repr__(self):
        return str(self)
    def T(self):
        return Matrix(zip(*self.data))
    def det(self, row=0):
        """Determinant using Laplace expansion on given row (slow)"""
        if len(self) == 1:
            return self[0,0]
        return sum((-1)**(i+row) * self[row,i]*self[[j for j in range(self.len[1])
                                                     if j != row],
                                                    [j for j in range(self.len[1])
                                                     if j != i]].det()
                   for i in range(len(self)))
    def map(self, f):
        return Matrix([[x if x is None else f(x) for x in y] for y in self.data], self.depth)
    def substitute(self, index, new):
        return self.map(lambda x: x if x is None else x.substitute(index, new))
    def satisfy(self, restrictions):
        """Returns a double list [[a,b],[c,d]...] where (a and b) or (c and d) or ... is
necessary/sufficient for given restrictions."""
        if restrictions is True:
            return [[self]]
        out = []
        for restriction in restrictions:
            to = [[self]]
            for r in restriction[3:5]:
                if not r is True:
                    to = [sum(and_, []) for x in to
                          for and_ in it.product(*[y.satisfy(r) for y in x])]
            if restriction[0] < self.depth or restriction[1] == '>':
                a = [[y.substitute(restriction[0],
                                   restriction[2] * (MRatio(MPoly({((0,)*restriction[0] + (1,)):1})) + 1))
                      for y in x]
                     for x in to]
            else:
                a = [[]]*len(to)
            if restriction[0] < self.depth or restriction[1] == '<':
                b = [[y.substitute(restriction[0],
                                   restriction[2] / (MRatio(MPoly({((0,)*restriction[0] + (1,)):1})) + 1))
                      for y in x]
                     for x in to]
            else:
                b = [[]]*len(to)
            out += [ap+bp for ap,bp in zip(a,b)]
        return out
    def solve_step(self, x, y):
        if self[x,y] is not None:
            raise ValueError("Already solved that index.")
        cpy = self[:,:]
        cpy.data[x][y] = MRatio(MPoly({((0,) * self.total + (1,)):1}))
        cpy.total += 1
        out = [[cpy]]
        for size in range(2, min(self.len) + 1):
            for rows in it.combinations(set(range(self.len[0])) - {x}, size - 1):
                truerows = sorted(rows + (x,))
                newrow = truerows.index(x)
                for cols in it.combinations(set(range(self.len[0])) - {y}, size - 1):
                    truecols = sorted(cols + (y,))
                    if all(x is not None for y in cpy[truerows,truecols] for x in y):
                        newout = []
                        for c in out:
                            newc = []
                            for d in c:
                                a = d[truerows,truecols].det(newrow).num.solve()
                                if a[0] is True:
                                    continue
                                if a[0] is False:
                                    break
                                if INFO:
                                    print(truerows, truecols)
                                    print(d)
                                    print(d[truerows,truecols].det(newrow))
                                    print(tuple(x for x in a[0] if x[0] == self.total))
                                    print()
                                newc.append(d.satisfy(tuple(x for x in a[0] if x[0] == self.total)))
                            else:
                                newout += [sum(x, []) for x in it.product(*newc)]
                        out = newout
        return out
INFO = False #print more information
a = MRatio({(1,):1})
b = MRatio({(0,1):1})
c = MRatio({(0,0,1):1})

A = Matrix([[1, 1, 1],
            [1, a+1, None],
            [1, (a+1) * (b+1), (a+1) * (b+1) * (c+1)]])
print('A', A.solve_step(1,2))

B = Matrix([[1, 1, 1],
            [1, None, a+1],
            [1, b+1, None]])
print('B', B.solve_step(1,1))

C = Matrix([[1, None, 1],
            [None, 1, 1],
            [1, a+1, (a+1) * (b+1)]])
print('C', C.solve_step(1,0))

D = Matrix([[1, None, 1],
            [None, None, 1],
            [1, 1, None]])
print('D', D.solve_step(1,0)[0][0].solve_step(0,1)[0][0].solve_step(1,1))
