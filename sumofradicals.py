from math import isqrt, gcd, lcm, prod, sumprod
from collections import defaultdict
from fractions import Fraction
from operator import mul
from sympy.ntheory import primefactors
from sympy.ntheory.factor_ import factorint
from itertools import product, islice, repeat
from functools import reduce, cache
from random import randint

import numpy as np
from scipy.linalg import circulant
from itertools import chain
from linalg import minor_laplace, det_laplace


@cache
def simplify_radical(n, r):
    """Return `s, n', r'` such that `sqrt[n](r)=s*sqrt[n'](r')`.
    
    The result has the lowest degree and radicand possible.
    `n` & `r` must be integers greater or equal 1.
    """
    #pull out perfect powers
    f = factorint(r)
    s = prod(pi**(ni//n) for pi, ni in f.items())
    f = {pi:ni%n for pi, ni in f.items()}
    #reduce degree
    while (cd := gcd(*f.values(), n)) > 1:
        f = {pi:ni//cd for pi, ni in f.items()}
        n //= cd
    r = prod(pi**ni for pi, ni in f.items())
    return s, n, r

#@total_ordering
class SumOfRadicals:
    r"""Class to handle numbers in the form $\sum_iv_i\sqrt[n_i]{r_i}/d$.
    
    The degrees and radicands are integers greater or equal 1,
    the factors in front of them are integers,
    the denominator is a integer greater or equal 1.
    This class is immutable.
    """
    
    #creation
    def __init__(self, n={}, d=1):
        """Create a new `SumOfRadicals`.
        
        The numerator `n` is expected to be a dictionary with `(ni, ri):vi` as
        keys:values, where the degrees 'ni' & radicands `ri`
        must be integers >=1 and the factors `vi` must be integers.
        `n` can also be an integer for an integer fraction.
        The denominator `d` must be a non-zero integer.
        No arguments defaults to 0.
        Numerator and denominator will be simplified to shortest terms.
        """
        if isinstance(n, int):
            n = {(1, 1):n}
        elif isinstance(n, dict):
            if not (all(isinstance(k, tuple) and len(k)==2 and isinstance(k[0], int) and isinstance(k[1], int) for k in n.keys())
                    and all(isinstance(v, int) for v in n.values())):
                raise TypeError('Degrees, radicands and factors must be integers.')
            if not all((n>=1 and r>=1) for n, r in n.keys()):
                raise ValueError('Degrees & Radicands must be greater than zero.')
        else:
            raise TypeError('Numerator must be an integer or a dictionary.')
        
        if isinstance(d, int):
            if not d:
                raise ValueError('Denominator must be non-zero.')
        else:
            raise TypeError('Denominator must be an integer.')
        
        #simplify numerator
        _n = defaultdict(int)
        for (n, r), v in n.items():
            s, n, r = simplify_radical(n, r)
            _n[(n, r)] += s * v
        n = {k:v for k, v in _n.items() if v}
        #short fraction
        while (cd := gcd(*n.values(), d)) > 1:
            n = {k:v//cd for k, v in n.items()}
            d //= cd
        #aesthetics
        if d < 0:
            n, d = {k:-v for k, v in n.items()}, -d
        n = dict(sorted(n.items()))
        
        self.n, self.d = n, d
    
    @staticmethod
    def random(N=5, precision=20):
        r"""Return a random `SumOfRadicals`.
        
        The factors from $\sqrt[1]{1}$ up to $\sqrt[N]{N}$ (incl.)
        and the denominator will be initialised with random integers `vi`
        such that `-precision//2 <= vi <= +precision//2`.
        """
        n = defaultdict(int, {(1, 1):randint(-precision//2, +precision//2)})
        for ni in range(2, N+1):
            for ri in range(2, N+1):
                n[(ni, ri)] = randint(-precision//2, +precision//2)
        #https://stackoverflow.com/a/69425586/7367030
        d = randint(-precision//2, +precision//2-1) or +precision//2
        return SumOfRadicals(n, d)
    
    
    #evaluation
    def __float__(self):
        """Calculate the `float` approximation.
        
        The `float` value might overflow.
        """
        #dicts aren't hashable, so __float__ can't be made a cached_property
        if not hasattr(self, '_float'):
            self._float = sumprod(self.values(), (r**(1/n) for n, r in self.keys())) \
                    / self.d
            #empty sumprod=int(0) nevertheless becomes float due to division
        return self._float
    
    def is_integer(self):
        """Return if this is an integer.
        
        If `True`, `int()` will return an integer,
        if `False`, `int()` will raise a `ValueError`.
        """
        return set(self.keys())<={(1, 1)} and self.d==1
    
    def __int__(self):
        """Return the integer value.
        
        If this doesn't represent an integer value a ValueError is raised.
        Can be checked beforehand with `is_integer()`.
        """
        if self.is_integer():
            return self.n.get((1, 1), 0)
        else:
            raise ValueError('doesn\'t represent an integer')
    
    def is_fraction(self):
        """Return this is an integer fraction.
        
        If `True`, `as_fraction()` will return a fraction,
        if `False`, `as_fraction()` will raise a `ValueError`.
        """
        return set(self.keys()) <= {(1, 1)}
    
    def as_fraction(self):
        """Return if this as an integer fraction.
        
        If this doesn't represent an integer fraction a ValueError is raised.
        Can be checked beforehand with `is_fraction()`.
        """
        if self.is_fraction():
            return Fraction(self.n.get((1, 1), 0), self.d)
        else:
            raise ValueError('doesn\'t represent a fraction')
    
    def __bool__(self):
        """Return if this is non-zero."""
        #leave this to avoid `if obj` call len(obj) instead of int(obj)
        #https://docs.python.org/3/reference/datamodel.html#object.__bool__
        #this is a numeric class, not a container
        try:
            return bool(int(self))
        except:
            return True
    
    
    #collection
    def __len__(self):
        """Return the number of summands in the numerator."""
        return len(self.n)
    
    def keys(self):
        """Return the radicands `r_i`."""
        return self.n.keys()
    
    def values(self):
        """Return the factors `n_i` infront of the square roots."""
        return self.n.values()
    
    def items(self):
        """Return the radicands `r_i` with factors `n_i` as tuples."""
        return self.n.items()
    
    
    #ordering
    def __eq__(self, other):
        """Return if this equals another `SumOfRadicals`, `int` or `Fraction`."""
        if isinstance(other, SumOfRadicals):
            return self.n == other.n and self.d == other.d
        elif isinstance(other, int) or isinstance(other, Fraction):
            try: #Fraction handles int comparison
                return self.as_fraction() == other
            except:
                return False
        else:
            raise TypeError
    
    def __abs__(self):
        """Return the absolute value as a `SumOfRadicals`."""
        return self if self>=0 else -self
    
    def __lt__(self, other):
        raise NotImplementedError
    
    
    #printing
    def __repr__(self):
        """Return a Unicode representation."""
        def int_to_superscript(n):
            return ''.join('⁰¹²³⁴⁵⁶⁷⁸⁹'[int(d)] for d in str(n))
        n = [f'{v:+d}{int_to_superscript(n)}{chr(0x221A)}{r}' for (n, r), v in self.items()]
        if len(n) <= 1: #no parentheses needed
            return (n[0] if n else '0') + f'/{self.d}'
        else:
            return '('+''.join(n)+')' + f'/{self.d}'
    
    def _repr_latex_(self):
        """Return a Latex representation."""
        n = [f'{v:+d}\\sqrt[{n}]{{{r}}}' for (n, r), v in self.items()]
        return '$\\frac{' + (''.join(n) if n else '0') + f'}}{{{self.d}}}$'
    
    
    #arithmetic
    #Implement the main arithmetic operations (+, *) with type checks
    #Inverse arithmetic (-, /) is completely reused
    def __add__(self, other):
        """Return the sum with an other `SumOfRadicals` or `int`."""
        if isinstance(other, SumOfRadicals):
            n = defaultdict(int)
            for k, v in self.items():
                n[k] += v * other.d
            for k, v in other.items():
                n[k] += v * self.d
            return SumOfRadicals(n, self.d*other.d)
        elif isinstance(other, int):
            n = defaultdict(int, self.n)
            n[(1, 1)] += other * self.d
            return SumOfRadicals(n, self.d)
        else:
            raise TypeError('can only add SumOfRadicals or int'
                    + f' (not "{type(other).__name__}") to SumOfRadicals')
    __radd__ = __add__
    
    def __neg__(self):
        """Return the additive negation."""
        return SumOfRadicals(self.n, -self.d)
    
    def __sub__(self, other):
        """Return the difference with an other `SumOfRadicals` or `int`."""
        return self + (-other)
    
    def __rsub__(self, other):
        """Return the difference from an other `SumOfRadicals` or `int`."""
        return (-self) + other
    
    def __mul__(self, other):
        """Return the product with an other `SumOfRadicals` or `int`."""
        if isinstance(other, SumOfRadicals):
            n = defaultdict(int)
            for (ni, ri), vi in self.items():
                for (nj, rj), vj in other.items():
                    nk = lcm(ni, nj)
                    r = ri**(nk//ni) * rj**(nk//nj)
                    n[(nk, r)] += vi * vj
            return SumOfRadicals(n, self.d*other.d)
        elif isinstance(other, int):
            return SumOfRadicals({k:v*other for k, v in self.items()}, self.d)
        else:
            raise TypeError('can only multiply SumOfRadicals or int'
                    + f' (not "{type(other).__name__}") with SumOfRadicals')
    __rmul__ = __mul__
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    @staticmethod
    def circulant_n(c, n):
        C = circulant(c)
        return n*np.triu(C, 1) + np.tril(C)
    
    @staticmethod
    def factorise_radical(n, r):
        """Return the prime factorisation of r with fractions (n as denominator) as exponents."""
        return {p:Fraction(e, n) for p, e in factorint(r).items()}
    
    @staticmethod
    def vectorise(s, degree, radical):
        c = [0]*degree
        for (n, r), v in s.items():
            f = SumOfRadicals.factorise_radical(n, r)
            if radical in f and degree%f[radical].denominator==0:
                i = f[radical].numerator * degree // f[radical].denominator
                c[i] += SumOfRadicals({(n, r//radical**(f[radical].numerator*n//f[radical].denominator)):v}, s.d)
            else:
                c[0] += SumOfRadicals({(n, r):v}, s.d)
        return c
    
    def __invert__(self):
        """Return the multiplicative reciprocal."""
        numerator, denominator = self.d, self*self.d
        while not denominator.is_fraction():
            factorisations = [SumOfRadicals.factorise_radical(n, r) for n, r in denominator.keys()]
            r = max(chain(*(f.keys() for f in factorisations)))
            n = lcm(*(f[r].denominator for f in factorisations if r in f))
            
            c = SumOfRadicals.vectorise(denominator, n, r)
            C = SumOfRadicals.circulant_n(c, r)
            C[0,:] = [SumOfRadicals({(n, r**i):1}) for i in range(C.shape[1])]
            factor = det_laplace(C)
            numerator *= factor
            denominator *= factor
        
        denominator = denominator.as_fraction()
        return numerator * denominator.denominator / denominator.numerator
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def __truediv__(self, other):
        """Return the quotient with an other `SumOfRadicals` or `int`."""
        if other == 0:
            raise ZeroDivisionError
        if isinstance(other, SumOfRadicals):
            return self * ~other
        elif isinstance(other, int):
            return SumOfRadicals(self.n, self.d*other)
        else:
            raise TypeError('can only divide SumOfRadicals by'
                    + f' SumOfRadicals or int (not "{type(other).__name__}")')
    
    def __rtruediv__(self, other):
        """Return an other `SumOfRadicals` or `int` divided by this."""
        return (~self) * other
    
    def __pow__(self, other):
        """Return this `SumOfRadicals` raised to some integer power.
        
        The exponent may be negative.
        """
        #repeat does typecheck for int
        if other >= 0:
            return reduce(mul, repeat(self, other), SumOfRadicals(1))
        else:
            return (~self)**(-other)



if __name__ == '__main__':
    from math import isclose as iscl
    
    def isclose(a, *b, rel_tol=1e-09, abs_tol=0.0):
        return all(iscl(a, bi, rel_tol=rel_tol, abs_tol=abs_tol) for bi in b)
    
    
    #creation & evaluation
    for _ in range(100):
        n = randint(-100, +100)
        a = SumOfRadicals(n)
        assert a==n and int(a)==n and isclose(float(a), n)
        
        n, d = randint(-100, +100), randint(-100, +100-1) or +100
        a = SumOfRadicals(n, d)
        assert a.as_fraction()==Fraction(n, d) and isclose(float(a), n/d)
    
    
    #comparison
    #for _ in range(100):
    #    a = SumOfRadicals.random(N=3)
    #    assert isclose(float(abs(a)), abs(float(a)))
    
    #for _ in range(100):
    #    a, b = SumOfRadicals.random(N=3), SumOfRadicals.random(N=3)
    #    assert (a<b) == (float(a)<float(b))
    
    
    #add
    for _ in range(1000):
        a, b = SumOfRadicals.random(), SumOfRadicals.random()
        assert isclose(float(a)+float(b), float(a+b))
        
        a, b = SumOfRadicals.random(), randint(-20, +20)
        assert isclose(float(a)+b, float(a+b), float(b+a))
    
    #sub
    for _ in range(1000):
        a, b = SumOfRadicals.random(), SumOfRadicals.random()
        assert isclose(float(a)-float(b), float(a-b))
        
        a, b = SumOfRadicals.random(), randint(-20, +20)
        assert isclose(float(a)-b, float(a-b), -float(b-a))
    
    #mul
    for _ in range(1000):
        a, b = SumOfRadicals.random(), SumOfRadicals.random()
        assert isclose(float(a)*float(b), float(a*b))
        
        a, b = SumOfRadicals.random(4), randint(-20, +20)
        assert isclose(float(a)*b, float(a*b), float(b*a))
    
    #invert
    for _ in range(10):
        a = SumOfRadicals.random(3)
        assert a / a == 1
    
    #div
    for _ in range(10):
        a, b = SumOfRadicals.random(5), SumOfRadicals.random(3)
        assert isclose(float(a)/float(b), float(a/b), rel_tol=1e-5)
    
        a, b = SumOfRadicals.random(5), randint(-20, +20-1) or +20
        assert isclose(float(a)/float(b), float(a/b))
    
        a, b = SumOfRadicals.random(3), randint(-20, +20)
        assert isclose(float(b)/float(a), float(b/a), rel_tol=1e-5)
    
    #pow
    for _ in range(10):
        a = SumOfRadicals.random(5)
        n = randint(0, +3)
        assert isclose(float(a**n), float(a)**n)
