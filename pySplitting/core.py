# AUTOGENERATED! DO NOT EDIT! File to edit: 00_core.ipynb (unless otherwise specified).

__all__ = ['findpos', 'isstrictlyascending', 'alpha', 'CreateMuVectors', 'CreateLyndonIndices', 'CreateEquation',
           'CreateConditions']

# Cell

from sympy import *
from sympy import sympify
from itertools import *


# Mathematica counts from 1, python from zero!!!!!!
def findpos(lam):
    """This determines whether the index string is a strictly ascending sequence.

    Parameters
    ----------
    lam : an indexable object of things that are comparable. Usually an array of integers

    Returns
    -------
    bool : the position in the string where two ints are equal.
    """
    retval = -1
    i = 0
    while (i < len(lam)-1) and (retval == -1):
        if lam[i] == lam[i+1]:
            i = i + 1
        else:
            retval = i
    if retval == -1:
        retval = len(lam)-1
    return retval

# Cell
def isstrictlyascending(lam):
    """This determines whether the index string is a strictly ascending sequence

    Parameters
    ----------
    lam : an indexable object of things that are comparable. Usually an array of integers

    Returns
    -------
    bool : true if the sequence is strictly ascending, else false.
    """
    retval = True
    for i in range(0, len(lam)-1):
#        print(lam[i],lam[i+1])
        retval = retval and lam[i] < lam[i+1]
#    i = 0
#    while (i < len(lam)-1) and retval == True:
#        if StrictLessThan(lam[i], lam[i+1]):
#        if sympify(lam[i] < lam[i+1]):
#            i = i+1
#        else:
#            retval = False
    return retval

# Cell
def alpha(lam):
    """
    The alpha factor. This determines the real alpha coefficient as given by Thalhammer and elaborated by Blanes.

    Parameters
    ----------
    lam : A multiindex of integers.

    Returns
    -------
    Rational: alpha(lam)
    """
    retval = 1
    if (len(lam) > 1) and not isstrictlyascending(lam):
        pos = findpos(lam)
        retval = 1/factorial(pos+1)*alpha(lam[pos+1:])
    return retval

# Cell
def CreateMuVectors(p,k):
    """This function creates the set of possible $$\mu$$ vectors.

    Parameters
    ----------
    p : the considered order
    k : 1 <=k<=p

    Returns
    -------
    array : A list of multiindices
    """
    absmu = p-k
    retval = []
    for el in product(range(0,absmu+1),repeat=k):
        if sum(el)<=absmu:
            retval.append(el)
    return retval

# Cell
def CreateLyndonIndices(p,k):
    """This function creates the set of Lyndon indices.

    Parameters
    ----------
    p : the considered order
    k : 1 <=k<=p

    Returns
    -------
    array : A list of Lyndon multiindices
    """
    retval = []
    S = list(range(1,p+1))#create the alphabet
    w = [S[0] - 1]
    while len(w) > 0:
        w[-1]=w[-1]+1
        m = len(w)
        if m == k:
            str = []
            for it in w:
                str.append(S[it-1])
            retval.append(str)
        # repeat w to get a string of length n
        while len(w) < k:
            w.append(w[-m])
        # remove the last character as long as it is equal to the largest character in S
        while (len(w) > 0) and w[-1] == S[-1]:
            w.pop()

    retval2=[]
    for el in retval:
        if sum(el)<=p:
            retval2.append([num - 1 for num in el])# subtract 1 from the final result. FIXME: this could be optimized
    return retval2

# Cell
def CreateEquation(mu, bvec, cvec):
    """This function gives the coefficient for a particular multiindex mu.

    It is the real coefficient of a Len(mu) long product of iterated commutators [A, B]_(mu_k).

    Parameters
    ----------
    mu : a multiindex
    bvec : a string of symbols.
    cvec : a string of symbols(usually the partial sums of the avec).

    Returns
    -------
    expr : a symbolic expression, a polynomial in terms of the contents of bvec and cvec.
"""
    k = len(mu)
    retval1=1
    for l in range(k):
        retval1=retval1*1/sympify(sum(mu[l:])+k-l)
    #print(retval1)
    retval2=0
    for it in combinations_with_replacement([*range(len(bvec),0,-1)],r=k):
        retvalprod = 1
        for i in range(k):
            retvalprod = retvalprod*bvec[it[i]-1]*cvec[it[i]-1]**mu[i]
        retval2=retval2 + alpha(it)*retvalprod

    return retval1-retval2

# Cell

from itertools import accumulate

def CreateConditions(p,avec,bvec,indexgenerator=CreateLyndonIndices):
    """This creates the set of equations using by default the Lyndon Basis elements.

    Parameters
    ----------
    p : the considered order
    avec: The set of symbols to use for the first operator.
    bvec: The set of symbols to use for the second operator.
    indexgenerator: (optional) by default we use indexgenerator for the Lyndon indices. Using CreateMuVectors
                     the indices from the overcomplete Hall-Basis can be used.

    Returns
    -------
    array : An array of Equations that have to be satisfied to fulfill the requested order p.
    """
    cvec=[*accumulate(avec)]
    cvec[-1]=1
    retval = [Eq(sum(avec),1)]
    for k in range(1,p+1):
        vecs=indexgenerator(p,k)
        for mu in vecs:
            retval.append(Eq(CreateEquation(mu,bvec,cvec),0))
    return retval