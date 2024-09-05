# sqrtfractions

Python module to handle linear combinations of square roots (field extension of rational numbers over integer square roots) exactly.
A single immutable class, `SqrtFraction`, is provided, to create object of the form

$$
    \frac{\sum_in_i\sqrt{r_i}}{d} \qquad n_i\in\mathbb{Z}, \ r_i, d\in\mathbb{N}^+
$$

where the values are represented by Pythons built-in arbitrary size integers, meaning there is no theoretical limit in magnitude nor precision. It's meant to be the next step from Pythons `fraction/Fraction` towards the reals.

## Installation

```console
pip install git+https://github.com/goessl/sqrtfractions.git
```

## Usage

A `SqrtFraction` can be initialised in two ways:
- with the constructor `SqrtFraction(n={}, d=1)`. The numerator `n` can be given as an integer or as a dictionary of keys:values that correspond to the radicand:factor terms.
- by the random factory `SqrtFraction.random(N=10, precision=20)`.
```python
>>> from sqrtfractions import SqrtFraction
>>> SqrtFraction(5)
+5√1/1
>>> SqrtFraction(5, 2)
+5√1/2
>>> SqrtFraction({2:3, 5:-7}, 2)
(+3√2-7√5)/2
```

`SqrtFraction`s can be printed
- in Unicode by `__repr__`
- in Latex by `_repr_latex_`.

`SqrtFraction`s can be casted to
- `float`s,
- `int`s (might check first with `is_integer()`),
- `Fraction`s (might check first with `is_fraction()`) &
- `bool`s.
```python
>>> float(SqrtFraction(5, 2))
2.5
>>> int(SqrtFraction(5))
5
>>> SqrtFraction(5, 2).as_fraction()
Fraction(5, 2)
>>> int(SqrtFraction(5, 2))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "sqrtfractions.py", line 114, in __int__
    raise ValueError('doesn\'t represent an integer')
ValueError: doesn't represent an integer
```

`SqrtFractions` are totally ordered
```python
>>> s, t = SqrtFraction({2:3}, 4), SqrtFraction({5:-6}, 7)
>>> s
+3√2/4
>>> t
-6√5/7
>>> s == t
False
>>> s > t
True
>>> abs(t)
+6√5/7
```

Basic arithmetic operations are implemented:
- unary negation `-` and inversion (reciprocal value) `~`,
- addition `+` and subtraction `-`, multiplication `*` and division `/` with other `SqrtFraction`s and with `int`s &
- exponentiation `**` to positive and negative integers.
```python
>>> -s
-3√2/4
>>> ~s
+2√2/3
>>> 
>>> s+t
(+21√2-24√5)/28
>>> s-t
(+21√2+24√5)/28
>>> 
>>> s*t
-9√10/14
>>> s/t
-7√10/40
>>> s**-2
+8√1/9
```

For more precise descriptions of the methods please refer to the docstrings.

## Design choices

- Reducing a new `SqrtFraction` directly after initialization: Initially `SqrtFraction`s didn't simplify themselves. It had been thought that the intermediate simplification during multiple consecutive operations would hinder performance. Eager simplification actually improved speed of an application. Actual testing hasn't been done.

## Relevant links

- [Rationalisation (mathematics) - Wikipedia](https://en.wikipedia.org/wiki/Rationalisation_(mathematics))
- [algorithms - Determine sign of sum of square roots - Mathematics Stack Exchange](https://math.stackexchange.com/a/1076510)

## TODO

- [ ] Hashing
- [ ] Arithmetic with floats (cast self to float and then use float arithmetic, like `fractions/Fraction`)
- [x] `abs` & `sign` methods. Might be difficult:
  [Square-root sum problem - Wikipedia](https://en.wikipedia.org/wiki/Square-root_sum_problem)
  [Sum of radicals - Wikipedia](https://en.wikipedia.org/wiki/Sum_of_radicals#:~:text=The%20sum%20of%20radicals%20is,finite%20linear%20combination%20of%20radicals%3A&text=are%20real%20numbers.)

## License (MIT)

Copyright (c) 2024 Sebastian Gössl

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
