from math import sqrt, gcd, sumprod
from collections import defaultdict
from random import randint
from itertools import product, islice
from sympy.ntheory import primefactors
from functools import cache, total_ordering


@cache
def _reduce_sqrt(r):
    """Return `s, t` such that `sqrt(r)=s*sqrt(t)`.
    
    The result has the lowest radicand possible.
    `r` must be a non-negative integer.
    """
    s, i = 1, 2 #factor in front, test factor
    while i**2 <= r:
        while not r % i**2:
            r //= i**2
            s *= i
        i += 1
    return s, r



@total_ordering
class SqrtFraction:
    r"""Class to handle numbers in the form $\sum_in_i\sqrt{r_i}/d$.
    
    The radicands are integers greater or equal 1,
    the factors infront of them are integers,
    the denominator is a non-zero integer.
    This class is immutable.
    """
    
    def __init__(self, n={}, d=1):
        """Construct a new SqrtFraction.
        
        The numerator `n` is expected to be a dictionary with r_i:n_i as
        keys:values where the radicands must be integers >=1
        and the factors must be integers.
        `n` can also be an integer for a usual fraction.
        The denominator `d` must be a non-zero integer.
        No arguments defaults to 0.
        """
        if isinstance(n, int):
            self.n, self.d = {1:n}, d
        elif isinstance(n, dict):
            assert all(isinstance(k, int) and k>0 for k in n.keys())
            assert all(isinstance(v, int) for v in n.values())
            assert isinstance(d, int) and d
            self.n, self.d = n, d
        else:
            raise TypeError
    
    @staticmethod
    def random(N=10, precision=20):
        r"""Return a random SqrtFraction.
        
        The factors from $\sqrt{1}$ up to $\\sqrt{N}$ (incl.)
        and the denominator will be initialised with random integers `ni`
        such that `-precision//2 <= ni <= +precision//2`.
        """
        n = {n:randint(-precision//2, +precision//2) for n in range(1, N+1)}
        #https://stackoverflow.com/a/69425586/7367030
        d = randint(-precision//2, +precision//2-1) or +precision//2
        return SqrtFraction(n, d)
    
    
    def __float__(self):
        """Calculate the float approximation.
        
        The float value might overflow. Then a copy of this object
        is first reduced and then calculating the float value is tried again.
        If this also overflows an OverflowError is raised.
        """
        try:
            #cast to float in case of and empty dict (then sum would be int)
            return float(sumprod(self.n.values(), map(sqrt, self.n.keys()))) \
                    / self.d
        except OverflowError:
            #return float(self.reduce()) may recourse indefinitely
            self = self.reduce()
            return float(sumprod(self.n.values(), map(sqrt, self.n.keys()))) \
                    / self.d
    
    def __int__(self):
        """Return the integer value.
        
        If this object doesn't represent an integer value
        a ValueError is raised.
        """
        self = self.reduce()
        if set(self.n.keys())<={1} and self.d==1:
            return self.n.get(1, 0)
        else:
            raise ValueError('doesn\'t represent an integer')
    
    
    def __len__(self):
        """Return the number of square roots in the numerator."""
        return len(self.n)
    
    def __eq__(self, other):
        """Return if this equals another SqrtFraction or integer in value."""
        if isinstance(other, SqrtFraction):
            self, other = self.reduce(), other.reduce()
            return self.n == other.n and self.d == other.d
        elif isinstance(other, int):
            try:
                return int(self) == other
            except:
                return False
        else:
            raise TypeError
    
    def __lt__(self, other):
        """Return if this is less than another SqrtFraction or integer."""
        #https://math.stackexchange.com/a/1076510
        if isinstance(other, SqrtFraction) or isinstance(other, int):
            s = (self - other).reduce()
            
            while len(s) > 1:
                #find highest factor in square roots
                pivot = max(
                        max(primefactors(n), default=0) for n in s.n.keys())
                
                #move all terms with factor to the right (s<?0 <=> (s-r)<?-r)
                r = SqrtFraction(
                        {k:-v for k, v in s.n.items() if k%pivot==0}, s.d)
                s = SqrtFraction(
                        {k:v for k, v in s.n.items() if k%pivot!=0}, s.d)
                
                #determine signs of both side for squaring
                #(inequality sign may change)
                #https://math.stackexchange.com/a/2347212
                #the pivot can be factored out on the right side
                #so sign(x) now has to deal with one degree less for both sides
                ss, sr = s<0, r*SqrtFraction({pivot:1}, pivot)<0
                #square both sides, adjust inequality sign
                s, r = s*s, r*r
                if ss:
                    s = -s
                if sr:
                    r = -r
                #pivot is now gone from within the square roots everywhere
                
                #move right side back to the left ((s-r)<?-r <=> s<?0)
                s = (s - r).reduce()
            
            #when there is only one term left in the numerator
            #then the sign is trivial
            return next(iter(s.n.values()), 0) * s.d < 0
        else:
            raise TypeError
    
    
    def __repr__(self):
        """Return a Unicode representation."""
        n = [f'{v:+d}{chr(0x221A)}{k}' for k, v in self.n.items()]
        return 'SqrtFraction{(' + ''.join(n) + f')/{self.d}}}'
    
    def _repr_latex_(self):
        """Return a Latex representation."""
        n = [f'{v:+d}\\sqrt{{{k}}}' for k, v in self.n.items()]
        return '$\\frac{' + (''.join(n) if n else '0') + f'}}{{{self.d}}}$'
    
    
    def reduce(self):
        """Return a mathematically simplified but equivalent `SqrtFraction`.
        
        The radicands in the square roots are reduced
        by pulling square factors out infront.
        Terms with a zero coefficient are removed.
        The fraction is shortened.
        A negative sign of the denominator is brought to the numerator.
        The square roots are sorted in ascending order.
        """
        #simplify numerator
        n = defaultdict(int)
        for k, v in self.n.items():
            f, k = _reduce_sqrt(k)
            n[k] += f * v
        n = {k:v for k, v in n.items() if v}
        #short fraction
        d = self.d
        while (cd := gcd(*n.values(), d)) > 1:
            n = {k:v//cd for k, v in n.items()}
            d //= cd
        #aesthetics
        if d < 0:
            n, d = {k:-v for k, v in n.items()}, -d
        return SqrtFraction(dict(sorted(n.items())), d)
    
    
    #Implement the main arithmetic operations (+, *) completely
    #with type checks and own integer calculations (don't just wrap the integer
    #into a SqrtFraction, this way it might be a little bit faster).
    #Just type checks for inverse arithmetic (-, /) to provide better
    #error strings.
    #Reuse commutative versions (radd, rmul) completely.
    def __add__(self, other):
        """Return the sum with an other `SqrtFraction` or `int`."""
        if isinstance(other, SqrtFraction):
            n = defaultdict(int)
            for k, v in self.n.items():
                n[k] += v * other.d
            for k, v in other.n.items():
                n[k] += v * self.d
            return SqrtFraction(n, self.d*other.d)
        elif isinstance(other, int):
            n = defaultdict(int, self.n)
            n[1] += other * self.d
            return SqrtFraction(n, self.d)
        else:
            raise TypeError('can only add SqrtFraction or int'
                    + f' (not "{type(other).__name__}") to SqrtFraction')
    __radd__ = __add__
    
    def __neg__(self):
        """Return the negation."""
        return SqrtFraction(self.n, -self.d)
    
    def __sub__(self, other):
        """Return the difference of an other `SqrtFraction` or `int`."""
        if isinstance(other, SqrtFraction) or isinstance(other, int):
            return self + (-other)
        else:
            raise TypeError('can only subtract SqrtFraction or int'
                    + f' (not "{type(other).__name__}") from SqrtFraction')
    
    def __rsub__(self, other):
        """Return the difference from an other `SqrtFraction` or `int`."""
        return (-self) + other
    
    def __mul__(self, other):
        """Return the product with an other `SqrtFraction` or `int`."""
        if isinstance(other, SqrtFraction):
            n = defaultdict(int)
            for ki, vi in self.n.items():
                for kj, vj in other.n.items():
                    n[ki*kj] += vi * vj
            return SqrtFraction(n, self.d*other.d)
        elif isinstance(other, int):
            return SqrtFraction({k:v*other for k, v in self.n.items()}, self.d)
        else:
            raise TypeError('can only multiply SqrtFraction or int'
                    + f' (not "{type(other).__name__}") with SqrtFraction')
    __rmul__ = __mul__
    
    def __invert__(self):
        """Return the reciprocal."""
        #https://www.youtube.com/watch?v=SjP6Mer0aL8
        #https://en.wikipedia.org/wiki/Rationalisation_(mathematics)
        t, r = {k:v for k, v in self.n.items() if v}, SqrtFraction(1)
        for p in islice(product((+1, -1), repeat=len(t)), 1, None):
            r *= SqrtFraction({k:s*v for (k, v), s in zip(t.items(), p)})
        return self.d * r / int(r * SqrtFraction(t))
    
    def __truediv__(self, other):
        """Return the quotient with an other `SqrtFraction` or `int`."""
        if isinstance(other, SqrtFraction):
            return self * ~other
        elif isinstance(other, int):
            return SqrtFraction(self.n, self.d*other)
        else:
            raise TypeError('can only divide SqrtFraction by'
                    + f' SqrtFraction or int (not "{type(other).__name__}")')
    
    def __rtruediv__(self, other):
        """Return the quotient from an other `SqrtFraction` or `int`."""
        return ~self * other



if __name__ == '__main__':
    assert _reduce_sqrt(0) == (1, 0)
    assert _reduce_sqrt(1) == (1, 1)
    assert _reduce_sqrt(2) == (1, 2)
    assert _reduce_sqrt(3) == (1, 3)
    assert _reduce_sqrt(4) == (2, 1)
    assert _reduce_sqrt(5) == (1, 5)
    assert _reduce_sqrt(6) == (1, 6)
    assert _reduce_sqrt(7) == (1, 7)
    assert _reduce_sqrt(8) == (2, 2)
    assert _reduce_sqrt(9) == (3, 1)
    assert _reduce_sqrt(10) == (1, 10)
    assert _reduce_sqrt(11) == (1, 11)
    assert _reduce_sqrt(12) == (2, 3)
    assert _reduce_sqrt(13) == (1, 13)
    assert _reduce_sqrt(14) == (1, 14)
    assert _reduce_sqrt(15) == (1, 15)
    assert _reduce_sqrt(16) == (4, 1)
    assert _reduce_sqrt(17) == (1, 17)
    assert _reduce_sqrt(18) == (3, 2)
    assert _reduce_sqrt(19) == (1, 19)
    assert _reduce_sqrt(20) == (2, 5)
    
    
    
    from math import isclose as iscl
    
    def isclose(a, *b, rel_tol=1e-09, abs_tol=0.0):
        return all(iscl(a, bi, rel_tol=rel_tol, abs_tol=abs_tol) for bi in b)
    
    
    #comparison
    for _ in range(1000):
        a, b = SqrtFraction.random(), SqrtFraction.random()
        assert (a<b) == (float(a)<float(b))
    
    #reduce
    for _ in range(1000):
        a = SqrtFraction.random()
        assert isclose(float(a), float(a.reduce()))
    
    #add
    for _ in range(1000):
        a, b = SqrtFraction.random(), SqrtFraction.random()
        assert isclose(float(a)+float(b), float(a+b))
        
        a, b = SqrtFraction.random(), randint(-20, +20)
        assert isclose(float(a)+b, float(a+b), float(b+a))
    
    #sub
    for _ in range(1000):
        a, b = SqrtFraction.random(), SqrtFraction.random()
        assert isclose(float(a)-float(b), float(a-b))
        
        a, b = SqrtFraction.random(), randint(-20, +20)
        assert isclose(float(a)-b, float(a-b), -float(b-a))
    
    #mul
    for _ in range(1000):
        a, b = SqrtFraction.random(), SqrtFraction.random()
        assert isclose(float(a)*float(b), float(a*b))
        
        a, b = SqrtFraction.random(4), randint(-20, +20)
        assert isclose(float(a)*b, float(a*b), float(b*a))
    
    #invert
    for _ in range(10):
        a = SqrtFraction.random(5)
        assert isclose(1/float(a), float(~a), rel_tol=1e-5)
    
    #div
    for _ in range(10):
        a, b = SqrtFraction.random(5), SqrtFraction.random(5)
        assert isclose(float(a)/float(b), float(a/b), rel_tol=1e-5)
        
        a, b = SqrtFraction.random(5), randint(-20, +20-1) or +20
        assert isclose(float(a)/float(b), float(a/b))
        
        a, b = SqrtFraction.random(5), randint(-20, +20)
        assert isclose(float(b)/float(a), float(b/a), rel_tol=1e-5)
