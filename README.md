# pySplitting
> A project to have python code that generates the order conditions for exponential splitting methods with the techniques of M. Thalhammer


This file will become your README and also the index of your documentation.

## Install

`pip pySplitting`

## How to use

Fill me in please! Don't forget code examples:

```python
alpha([1,1,1])
```




$\displaystyle \frac{1}{6}$


from sympy import symarray

avec=[*symarray('a',3)]
bvec=[*symarray('b',3)]
print(CreateConditions(3,avec,bvec))
print("=========================")
print(CreateConditions(3,avec,bvec,CreateMuVectors))