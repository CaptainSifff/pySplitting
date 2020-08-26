# pySplitting
> A project to have python code that generates the order conditions for exponential splitting methods with the techniques of M. Thalhammer.


## Install

`pip pySplitting`

## How to use

This is a python project that generates the equations for the order conditions such that exponential splitting methods of the given order or with the desired features can be designed. To that end we utilize sympy to get symbolic expressions that can be further processed.

First we will derive the equations for the familiar leapfrog/Strang splitting. This is a method of order p=2 and uses the indeterminates ((t0,t1),(v0,v1))=((1/2,1/2).(1,0)). Hence we have to set up sympy first to give us the respective symbols:

```python
from sympy.matrices.dense import symarray
from sympy import *
tvec=[*symarray('t',2)]
vvec=[*symarray('v',2)]
```

Next we build the equations. To that end we use the function CreateCondition().
The first argument is the order, the second and third argument are the indeterminates for the resulting polynomials.

```python
strangconds=CreateConditions(2,tvec,vvec)
display(strangconds)
```


    [Eq(t_0 + t_1, 1), Eq(-v_0 - v_1 + 1, 0), Eq(-t_0*v_0 - v_1 + 1/2, 0)]


Now we check wether the solutions fulfill these equations:

```python
[it.subs([('t_0',sympify(1)/2),('t_1',sympify(1)/2),('v_0',1),('v_1',0)]) for it in strangconds]
```




    [True, True, True]



It works! Now we move on to check the improved Leapfrog by McLachlan, and utilize that we can impose already some structure(symmetry and the consistency equations) on the indeterminates to get some simplifications:

```python
tvec=[*symarray('t',3)]
tvec[2]=tvec[0]
tvec[1]=1-2*tvec[0]
vvec=[*symarray('v',3)]
vvec[2]=0
vvec[1]=vvec[0]
print("Indeterminates: ",tvec,vvec)
McLachlanconds=CreateConditions(2,tvec,vvec)
print(McLachlanconds)
```

    Indeterminates:  [t_0, 1 - 2*t_0, t_0] [v_0, v_0, 0]
    [True, Eq(1 - 2*v_0, 0), Eq(-t_0*v_0 - v_0*(1 - t_0) + 1/2, 0)]


We check that the given numbers(v_0=1/2,t_0=0.1931) satisfy the equation:

```python
[it.subs([('t_0',0.19318332),('v_0',sympify(1)/2)]) for it in McLachlanconds]
```




    [True, True, True]



True! Now we move on to two more fourth order esamples from the literature:
