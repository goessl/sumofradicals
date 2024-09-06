# sumofradicals

Python module to handle linear combinations of radicals (field extension of rational numbers over integer roots of various degree) exactly.
A single immutable class, `SumOfRadicals`, is provided, to create object of the form

$$
    \frac{\sum_iv_i\sqrt[n_i]{r_i}}{d} \qquad v_i\in\mathbb{Z}, \ n_i, r_i, d\in\mathbb{N}^+
$$

where the values are represented by Pythons built-in arbitrary size integers, meaning there is no theoretical limit in magnitude nor precision. It's meant to be the next step from Pythons `fraction/Fraction` towards the reals.

*A clearer and faster version for square roots only is available as [sqrtfractions](https://github.com/goessl/sqrtfractions).*

## Installation

```console
pip install git+https://github.com/goessl/sumofradicals.git
```

## Usage

A `SumOfRadicals` can be initialised in two ways:
- with the constructor `SumOfRadicals(n={}, d=1)`. The numerator `n` can be given as an integer or as a dictionary of keys:values that correspond to the (degree, radicand):factor terms.
- by the random factory `SumOfRadicals.random(N=10, precision=20)`.
```python
>>> SumOfRadicals(5)
+5¹√1/1
>>> SumOfRadicals(5, 2)
+5¹√1/2
>>> SumOfRadicals({(1, 1):1, (2, 3):11, (5, 7):13}, 17)
(+1¹√1+11²√3+13⁵√7)/17
```

`SumOfRadicals`s can be printed
- in Unicode by `__repr__`
- in Latex by `_repr_latex_`.

`SumOfRadicals`s can be casted to
- `float`s,
- `int`s (might check first with `is_integer()`),
- `Fraction`s (might check first with `is_fraction()`) &
- `bool`s.
```python
>>> float(SumOfRadicals(5, 2))
2.5
>>> int(SumOfRadicals(5))
5
>>> SumOfRadicals(5, 2).as_fraction()
Fraction(5, 2)
>>> int(SumOfRadicals(5, 2))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "sumofradicals.py", line 134, in __int__
    raise ValueError('doesn\'t represent an integer')
ValueError: doesn't represent an integer
```

Basic arithmetic operations are implemented:
- unary negation `-`,
- addition `+` and subtraction `-` and multiplication `*` with other `SumOfRadicals`s and with `int`s &
- exponentiation `**` to non-negative integers.
```python
>>> s, t = SumOfRadicals({(2, 3):5}, 7), SumOfRadicals({(1, 1):2}, 11)
>>> s
+5²√3/7
>>> t
+2¹√1/11
>>> -s
-5²√3/7
>>> s+t
(+14¹√1+55²√3)/77
>>> s-t
(-14¹√1+55²√3)/77
>>> s*t
+10²√3/77
>>> s**2
+75¹√1/4
```

For more precise descriptions of the methods please refer to the docstrings.

## Design choices

- Fractions instead of a purely integer linear combinations: In this form all basic arithmetic operations are closed (division wouldn't be otherwise).

## TODO

- [ ] Hashing
- [ ] Arithmetic with floats (cast self to float and then use float arithmetic, like `fractions/Fraction`)
- [ ] `abs` & `sign` methods. Might be difficult:
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
