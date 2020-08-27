# pySplitting
> A project to have python code that generates the order conditions for exponential splitting methods with the equations of M. Thalhammer given in the 2008 paper "High-Order Exponential Operator Splitting Methods for Time-Dependent Schrödinger Equations" available at ´https://epubs.siam.org/doi/abs/10.1137/060674636?mobileUi=0´


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



It works! Now we move on to check the improved Leapfrog by McLachlan(McLachlan 1995, https://www.massey.ac.nz/~rmclachl/sisc95.pdf), and utilize that we can impose already some structure(symmetry and the consistency equations) on the indeterminates to get some simplifications:

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
First the fourth order method from McLachlan's 1995 paper:

```python
tvec=[*symarray('t',5)]
tvec[4]=tvec[0]
tvec[3]=tvec[1]
tvec[2]=1-2*(tvec[0]+tvec[1])
vvec=[*symarray('v',5)]
vvec[4]=0
vvec[3]=vvec[0]
vvec[1]=(1-2*vvec[0])/2
vvec[2]=vvec[1]
print("Indeterminates: ",tvec,vvec)
McLachlanconds4=CreateConditions(4,tvec,vvec)
print(McLachlanconds4)
[it.subs([('v_0',S(6)/S(11)),('t_0',(642+Pow(471,S(1)/2))/S(3924)),('t_1',S(121)/S(3924)*(12-Pow(471,S(1)/2)))]) for it in McLachlanconds4]
```

    Indeterminates:  [t_0, t_1, -2*t_0 - 2*t_1 + 1, t_1, t_0] [v_0, 1/2 - v_0, 1/2 - v_0, v_0, 0]
    [True, True, Eq(-t_0*v_0 - v_0*(1 - t_0) - (1/2 - v_0)*(t_0 + t_1) - (1/2 - v_0)*(-t_0 - t_1 + 1) + 1/2, 0), Eq(-t_0**2*v_0 - v_0*(1 - t_0)**2 - (1/2 - v_0)*(t_0 + t_1)**2 - (1/2 - v_0)*(-t_0 - t_1 + 1)**2 + 1/3, 0), Eq(-t_0**3*v_0 - v_0*(1 - t_0)**3 - (1/2 - v_0)*(t_0 + t_1)**3 - (1/2 - v_0)*(-t_0 - t_1 + 1)**3 + 1/4, 0), Eq(-3*t_0*v_0**2/2 - 2*t_0*v_0*(1/2 - v_0) - v_0**2*(1 - t_0)/2 - v_0*(1/2 - v_0)*(t_0 + t_1) - v_0*(1/2 - v_0)*(-t_0 - t_1 + 1) - 3*(1/2 - v_0)**2*(t_0 + t_1)/2 - (1/2 - v_0)**2*(-t_0 - t_1 + 1)/2 + 1/6, 0), Eq(-3*t_0**2*v_0**2/2 - 2*t_0**2*v_0*(1/2 - v_0) - v_0**2*(1 - t_0)**2/2 - v_0*(1/2 - v_0)*(t_0 + t_1)**2 - v_0*(1/2 - v_0)*(-t_0 - t_1 + 1)**2 - 3*(1/2 - v_0)**2*(t_0 + t_1)**2/2 - (1/2 - v_0)**2*(-t_0 - t_1 + 1)**2/2 + 1/12, 0), Eq(-7*t_0*v_0**3/6 - 3*t_0*v_0**2*(1/2 - v_0) - 2*t_0*v_0*(1/2 - v_0)**2 - v_0**3*(1 - t_0)/6 - v_0**2*(1/2 - v_0)*(t_0 + t_1)/2 - v_0**2*(1/2 - v_0)*(-t_0 - t_1 + 1)/2 - 3*v_0*(1/2 - v_0)**2*(t_0 + t_1)/2 - v_0*(1/2 - v_0)**2*(-t_0 - t_1 + 1)/2 - 7*(1/2 - v_0)**3*(t_0 + t_1)/6 - (1/2 - v_0)**3*(-t_0 - t_1 + 1)/6 + 1/24, 0)]





    [True, True, True, True, True, True, True, True]



Next we check the popular S_6 4 method by from Blanes2002("Practical symplectic partitioned Runge–Kutta and Runge–Kutta–Nyström methods")
Since it's a symmetric method it is sufficient to work with the third order equations.

```python
tvec=[*symarray('t',7)]
tvec[6]=tvec[0]
tvec[5]=tvec[1]
tvec[4]=tvec[2]
tvec[3]=1-2*(tvec[0]+tvec[1]+tvec[2])
vvec=[*symarray('v',7)]
vvec[2]=(1-2*(vvec[0]+vvec[1]))/2
vvec[6]=0
vvec[5]=vvec[0]
vvec[4]=vvec[1]
vvec[3]=vvec[2]
print("Indeterminates: ",tvec,vvec)
S64conds=CreateConditions(3,tvec,vvec)
display(S64conds)
```

    Indeterminates:  [t_0, t_1, t_2, -2*t_0 - 2*t_1 - 2*t_2 + 1, t_2, t_1, t_0] [v_0, v_1, -v_0 - v_1 + 1/2, -v_0 - v_1 + 1/2, v_1, v_0, 0]



    [True,
     True,
     Eq(-t_0*v_0 - v_0*(1 - t_0) - v_1*(t_0 + t_1) - v_1*(-t_0 - t_1 + 1) - (t_0 + t_1 + t_2)*(-v_0 - v_1 + 1/2) - (-v_0 - v_1 + 1/2)*(-t_0 - t_1 - t_2 + 1) + 1/2, 0),
     Eq(-t_0**2*v_0 - v_0*(1 - t_0)**2 - v_1*(t_0 + t_1)**2 - v_1*(-t_0 - t_1 + 1)**2 - (t_0 + t_1 + t_2)**2*(-v_0 - v_1 + 1/2) - (-v_0 - v_1 + 1/2)*(-t_0 - t_1 - t_2 + 1)**2 + 1/3, 0),
     Eq(-3*t_0*v_0**2/2 - 2*t_0*v_0*v_1 - 2*t_0*v_0*(-v_0 - v_1 + 1/2) - v_0**2*(1 - t_0)/2 - v_0*v_1*(t_0 + t_1) - v_0*v_1*(-t_0 - t_1 + 1) - v_0*(t_0 + t_1 + t_2)*(-v_0 - v_1 + 1/2) - v_0*(-v_0 - v_1 + 1/2)*(-t_0 - t_1 - t_2 + 1) - 3*v_1**2*(t_0 + t_1)/2 - v_1**2*(-t_0 - t_1 + 1)/2 - 2*v_1*(t_0 + t_1)*(-v_0 - v_1 + 1/2) - v_1*(t_0 + t_1 + t_2)*(-v_0 - v_1 + 1/2) - v_1*(-v_0 - v_1 + 1/2)*(-t_0 - t_1 - t_2 + 1) - 3*(t_0 + t_1 + t_2)*(-v_0 - v_1 + 1/2)**2/2 - (-v_0 - v_1 + 1/2)**2*(-t_0 - t_1 - t_2 + 1)/2 + 1/6, 0)]


And now let's see whether it works

```python
res=[it.simplify().subs([
    ('t_0',S(0.0792036964311956)),
          ('t_1',S(0.35317290604977410)),
          ('t_2',S(-0.0420650803577195)),
          ('v_0',S(0.2095151066133620)),
          ('v_1',S(-0.143851773179818))
         ]) for it in S64conds]
print(res)
```

    [True, True, True, False, False]


```python
print(S64conds)
```

    [True, True, Eq(-t_0*v_0 - v_0*(1 - t_0) - v_1*(t_0 + t_1) - v_1*(-t_0 - t_1 + 1) - (t_0 + t_1 + t_2)*(-v_0 - v_1 + 1/2) - (-v_0 - v_1 + 1/2)*(-t_0 - t_1 - t_2 + 1) + 1/2, 0), Eq(-t_0**2*v_0 - v_0*(1 - t_0)**2 - v_1*(t_0 + t_1)**2 - v_1*(-t_0 - t_1 + 1)**2 - (t_0 + t_1 + t_2)**2*(-v_0 - v_1 + 1/2) - (-v_0 - v_1 + 1/2)*(-t_0 - t_1 - t_2 + 1)**2 + 1/3, 0), Eq(-3*t_0*v_0**2/2 - 2*t_0*v_0*v_1 - 2*t_0*v_0*(-v_0 - v_1 + 1/2) - v_0**2*(1 - t_0)/2 - v_0*v_1*(t_0 + t_1) - v_0*v_1*(-t_0 - t_1 + 1) - v_0*(t_0 + t_1 + t_2)*(-v_0 - v_1 + 1/2) - v_0*(-v_0 - v_1 + 1/2)*(-t_0 - t_1 - t_2 + 1) - 3*v_1**2*(t_0 + t_1)/2 - v_1**2*(-t_0 - t_1 + 1)/2 - 2*v_1*(t_0 + t_1)*(-v_0 - v_1 + 1/2) - v_1*(t_0 + t_1 + t_2)*(-v_0 - v_1 + 1/2) - v_1*(-v_0 - v_1 + 1/2)*(-t_0 - t_1 - t_2 + 1) - 3*(t_0 + t_1 + t_2)*(-v_0 - v_1 + 1/2)**2/2 - (-v_0 - v_1 + 1/2)**2*(-t_0 - t_1 - t_2 + 1)/2 + 1/6, 0)]


## Complex Hermitian Methods where one set of coefficients is real

In this section we check our new hermitian methods where one set of coefficients is real.

### CHR_3 3

```python
avec=[*symarray('a',4)]
bvec=[*symarray('b',4)]
tvec=[
    avec[0]+I*bvec[0],
    avec[1]+I*bvec[1],
    avec[1]-I*bvec[1],
    avec[0]-I*bvec[0],
]
vvec=[*symarray('v',4)]
vvec[1]=1-2*vvec[0]
vvec[3]=0
vvec[2]=vvec[0]
print("Indeterminates: ",tvec,vvec)
CHR33conds=CreateConditions(3,tvec,vvec)
display(CHR33conds)
```

    Indeterminates:  [a_0 + I*b_0, a_1 + I*b_1, a_1 - I*b_1, a_0 - I*b_0] [v_0, 1 - 2*v_0, v_0, 0]



    [Eq(2*a_0 + 2*a_1, 1),
     True,
     Eq(-v_0*(a_0 + I*b_0) - v_0*(a_0 + 2*a_1 + I*b_0) - (1 - 2*v_0)*(a_0 + a_1 + I*b_0 + I*b_1) + 1/2, 0),
     Eq(-v_0*(a_0 + I*b_0)**2 - v_0*(a_0 + 2*a_1 + I*b_0)**2 - (1 - 2*v_0)*(a_0 + a_1 + I*b_0 + I*b_1)**2 + 1/3, 0),
     Eq(-3*v_0**2*(a_0 + I*b_0)/2 - v_0**2*(a_0 + 2*a_1 + I*b_0)/2 - v_0*(1 - 2*v_0)*(a_0 + I*b_0) - v_0*(1 - 2*v_0)*(a_0 + a_1 + I*b_0 + I*b_1) - (1 - 2*v_0)**2*(a_0 + a_1 + I*b_0 + I*b_1)/2 + 1/6, 0)]


```python
res=[it.simplify().subs([
    ('a_0',3*S(1)/S(24)),
          ('a_1',3*S(1)/S(8)),
      ('b_0',-Pow(3,S(1)/S(2))/S(24)),
          ('b_1',Pow(3,S(1)/S(2))/S(8)),
          ('v_0',S(1)/S(3))
         ]).simplify() for it in S64conds]
print(res)
```

    [True, True, True, True, True]


### CHR_5 4

```python
avec=[*symarray('a',6)]
bvec=[*symarray('b',6)]
tvec=[
    avec[0]+I*bvec[0],
    avec[1]+I*bvec[1],
    avec[1]-I*bvec[1],
    avec[0]-I*bvec[0],
]
print("Indeterminates: ",tvec,vvec)
CHR54conds=CreateConditions(4,tvec,vvec)
display(CHR54conds)
```

### CHR_9 5

```python
avec=[*symarray('a',10)]
bvec=[*symarray('b',10)]
tvec=[
    avec[0]+I*bvec[0],
    avec[1]+I*bvec[1],
    avec[1]-I*bvec[1],
    avec[0]-I*bvec[0],
]
print("Indeterminates: ",tvec,vvec)
CHR95conds=CreateConditions(5,tvec,vvec)
display(CHR95conds)
```

### CHR_15 6

```python
avec=[*symarray('a',16)]
bvec=[*symarray('b',16)]
tvec=[
    avec[0]+I*bvec[0],
    avec[1]+I*bvec[1],
    avec[1]-I*bvec[1],
    avec[0]-I*bvec[0],
]
print("Indeterminates: ",tvec,vvec)
CHR156conds=CreateConditions(6,tvec,vvec)
display(CHR156conds)
```
