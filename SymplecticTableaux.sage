# ****************************************************************************
#       Copyright (C) 2023 Álvaro Gutiérrez <gutierrez.caceres@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

r"""
Implementation of symplectic tableaux and symplectic patterns.

Symplectic tableaux are type C analogues of semistandard Young tableaux.
Explicitly, the generating function of symplectic tableaux of a given shape
`\lambda` (a partition) is the character of the irreducible representation
of `\mathfrak{sp}(2n)` indexed by `\lambda`. 

In the literature, several definitions of symplectic tableaux are available.
We implement four of these; King's tableaux [Kin]_, De Concini's tableaux [DeC]_,
Krattenthaler's tableaux [Kra]_, and Kashiwara's tableaux [NK]_. Note that some
functions on Kashiwara's tableaux are already implemented in Sage as part of the
"tensor product of crystals" library.


AUTHORS

-  Álvaro Gutiérrez (2023): initial version

.. rubric:: REFERENCES

.. [DeC79] Corrado De Concini.
    Symplectic standard tableaux.
    Advances in Mathematics, 34(1):1–27, (1979).
    :doi:`10.1016/0001-8708(79)90061-6`

.. [Kin76] Ronald King.
    Weight multiplicities for the classical groups.
    Group Theoretical Methods in Physics, 490–499, (1976).
    
.. [KN94] Masaki Kashiwara and Toshiki Nakashima.
    Crystal Graphs for Representations of the q-Analogue of Classical Lie Algebras.
    Journal of Algebra, 165(2):295–345, (1994). 
    :doi:`10.1006/jabr.1994.1114`
   
.. [Kra98] Christian Krattenthaler.
    Identities for classical group characters of nearly rectangular shape.
    Journal of Algebra, 209(1):1–64, (1998).
    :doi:`10.1006/jabr.1998.7531`
    
.. [She99] Jeffrey Sheats.
    A symplectic jeu de taquin bijection between the tableaux of King and of De Concini.
    Transactions of the American Mathematical Society, 351(9):3569–3607, (1999).
    http://www.jstor.org/stable/117933
"""

class SymplecticTableau(SageObject):
    r"""
    Four different combinatorial models for symplectic tableaux
    are implemented. Internally, these tableaux are always in the
    alphabet ``1 <\cdots< \cdots < 2n.`` Later, they may be displayed
    with a different alphabet. 
    
    The default model is King's. They are displayed in the alphabet
    ``1 < 1' <\cdots< 2' < \cdots <\cdots< n'.``
    They are well-defined if every column is admissible, that is,
    if for every `i` there are at most `i` numbers smaller or equal
    to `i'` in each column.
    
    A second model is De Concini's. These are displayed in the alphabet
    ``n' < \cdots < 2' < 1' <\cdots< 2 < \cdots < n.``
    These are well-defined if every column is admissible (after
    reordering to the King's alphabet) and their split version is
    semistandard. See [She99]
    
    A third model is Krattenthaler's. These are displayed in the
    alphabet ``1 <\cdots< \cdots < 2n.``
    These are well-defined when they are the split version of a
    De Concini tableau.
    
    A final model is Kashiwara's. These are displayed in the alphabet
    ``1 <\cdots< \cdots <\cdots< n' < \cdots < 2' < 1'.``
    These are well-defined when they are their cosplit version is
    a Krattenthaler tableau.
    
    Additionally, a custom type of tableaux can be defined by specifying
    the alphabet in which they should be displayed.
    
    Skew tableaux are implemented by introducing 0s.
    
    EXAMPLES::
    sage: tab = SymplecticTableau(3, rows = [[1,1,2,3,3],[4,4,5,6]]); tab
    1  1  1' 2' 2'
    2  2  3  3'
    sage: tab.is_well_defined()
    True
    sage: tab.DeConcini()
    3' 2' 2' 2  2
    1' 1' 1  3
    sage: tab.Kashiwara()
    1  2  2  2' 2'
    3  3  3' 1'
    sage: tab.f(1)
    1  1' 1' 2' 2'
    2  2  3  3'
    sage: print(tab)
    Symplectic tableau of type King, shape [5, 4], and weight x1^(-1)x2^(0)x3^(1)
    sage: tab.shape()
    [5, 4]
    sage: latex(tab)
    \ytableaushort{{1 }{1 }{1'}{2'}{2'},{2 }{2 }{3 }{3'}}
    """
    def __init__(self, n, rows = None, cols = None, type = 'King', alphabet = None):
        r"""
        INPUT:
            - n : Positive integer
            - rows : SemistandardTableau or list of lists of integers
            - cols : SemistandardTableau or list of lists of integers
            - type : Either 'King', 'DeConcini', 'Krattenthaler' or 'Kashiwara'. (Default: 'King')
            - alphabet : a list of strings. (Default: None)
        """
        assert bool(rows == None) != bool(cols == None), "Exactly one of 'rows' or 'cols' must be provided"
        assert bool(type == None) != bool(alphabet == None), "Exactly one of 'type' or 'alphabet' must be provided"
        
        if rows != None:
            self._rows = [list(row) for row in rows]
            self._cols = transpose(self._rows)
        else:
            self._cols = [list(col) for col in cols]
            self._rows = transpose(self._cols)
        self._t = SemistandardTableau(self._rows)
        self._len = len(self._rows)
        self._n = n
        if alphabet == None:
            self._type = type
            if type == 'King':
                self._alphabet = sum([[str(i), str(i) + "'"] for i in [1..n]], [])
            elif type == 'DeConcini':
                self._alphabet = [str(i) + "'" for i in [1..n][::-1]] + [str(i) for i in [1..n]]
            elif type == 'Krattenthaler' or type == "split":
                self._alphabet = [str(i) for i in [1..2*n]]
            elif type == 'Kashiwara':
                self._alphabet = [str(i) for i in [1..n]] + [str(i) + "'" for i in [1..n][::-1]]
        else:
            self._type = 'Custom'
            self._alphabet = alphabet
        self._shape = self.shape()
        self._weight = self.weight()
    def __eq__(self, other):
        return self.King()._rows == other.King()._rows
    def __repr__(self):
        if self._rows == []:
            return ''
        
        l = floor(log(self._n)/log(10))+2
        alph = ['*'] + self._alphabet
        alph2 = []
        for i in range(len(alph)):
            alph2.append(alph[i] + " "*(l - len(alph[i])))
        return '\n'.join(' '.join(alph2[i] for i in row) for row in self._rows)
    def __str__(self):
        if self._type == "Krattenthaler" or self._type == "split" or self._type == "Custom":
            return f"Symplectic tableau of type {self._type} and shape {self._shape}"
        return f"Symplectic tableau of type {self._type}, shape {self._shape}, and weight {self._weight}"
    def cols(self):
        return self._cols
    def transpose(self):
        return self._cols
    def __iter__(self):
        return SymplecticTableauIterator(self)
    def __len__(self):
        return self._len
    def _latex_(self):
        alph = [''] + self._alphabet
        return r"\ytableaushort{" \
                + ','.join('{' + '}{'.join(alph[i] for i in row) + '}' \
                for row in self._rows) \
                + '}'
    def list(self):
        return self._rows
    def rows(self):
        return self._rows
    def len(self):
        return self._len
    def dict(self):
        return {i: self._rows[i] for i in self._len}
    def is_well_defined(self):
        r"""
        Checks if the tableau ``self`` is a tableau of type ``self._type``.
        """
        if self._type == 'King':
            return all(self._rows[i][0] >= 2*i+1 for i in range(len(self._rows)))
        elif self._type == 'DeConcini':
            return all([is_admissible(c, self._n) for c in self._cols]) and splitVersion(self)._t.is_semistandard()
        elif self._type == 'Krattenthaler' or self._type == 'split':
            return cosplitVersion(cosplitInverse(self)) == self and self._t.is_semistandard()
        elif self._type == 'Kashiwara':
            return all([is_coadmissible(c, self._n) for c in self._cols]) and cosplitVersion(self)._t.is_semistandard()
        raise ValueError("Well-definedness cannot be checked on a custom tableau")
    def n(self):
        return self._n
    def size(self):
        return self._n
    def weight(self):
        if self._type == "Krattenthaler" or self._type == "split" or self._type == "Custom":
            return "Weight is only defined for King, De Cocini, or Kashiwara tableaux."
        l = len(self._alphabet[0])
        alph = ['*'] + self._alphabet
        R = sum([[alph[i] for i in row] for row in self._rows], [])
        
        w = {}
        for i in [1..self._n]:
            w[i] = R.count("%s"%i) - R.count("%s'"%i)
        return ''.join('x%s^(%s)'%(i, w[i]) for i in [1..self._n] if w[i] != 0)
    def shape(self):
        if self._type == 'Krattenthaler' or self._type == "split":
            return Partition([int(len(row)/2) for row in self._rows])
        else:
            return Partition([len(row) for row in self._rows])
    def bender_knuth_involution(self, i, display=False):
        r"""
        Performs the `i`th Bender--Knuth involution on the tableau.

        The `i`th Bender--Knuth involution applied to a tableau with
        weight `x^{\lambda}` returns a tableau of same shape and 
        with weight `s_i.x^{\lambda}`, where `s_0, ..., s_{n-1}` are the
        generators of the Weyl group of type C. These are analogs of type
        A Bender--Knuth involutions.
        
        To compute the `0`th Bender--Knuth involution on a King tableau,
        we interpret it as a semistandard tableau in the alphabet `1<\cdots<2n`
        and then perform the `1`st (type A) Bender--Knuth involution.
        
        For `i = 1, ..., n-1`, the `i`th Bender--Knuth involution on a King
        tableau is defined as the composite of five maps. The first four maps
        are usual type A Bender--Knuth involutions. More precisely, if the
        tableau is interpreted as a (type A) semistandard tableau in the alphabet
        `1<\cdots<2n`, then the first four maps are the `2i`th, the `(2i-1)`st, the
        `(2i+1)`st, and the `2i`th Bender--Knuth involutions. The fifth map 
        changes every `\{i, i'\}`-vertical domino between rows `i` and `i+1` for
        a `\{i+1, i+1'\}`-vertical domino. It then resorts rows `i` and `i+1` as
        to make them increasing.
        
        The function converts the tableau to a King tableau, performs the involution
        and converts back to the original type. If ``display`` is set to ``True``,
        then the intermediate stages of the algorithm are displayed.
        
        EXAMPLES::
            sage: T = SymplecticTableau(3, rows = [[1,3,4],[5,5,6],[6,6]]); T
            1  2  2'
            3  3  3'
            3' 3'
            sage: T.bender_knuth_involution(0)
            1' 2  2'
            3  3  3'
            3' 3'
            sage: T.bender_knuth_involution(2, display = True)
            Now applying Bender--Knuth to the tableau
            T0 :
            1  2  2'
            3  3  3'
            3' 3'
            ---------------
            T1 :
            1  2  3
            2' 2' 3'
            3' 3'
            ---------------
            T2 :
            1  2  3
            2  2' 3'
            3' 3'
            ---------------
            T3 :
            1  2  3
            2  2' 3'
            3  3
            ---------------
            T4 :
            1  2  2'
            2  2' 3'
            2' 3
            ---------------
            T5 :
            1  2  2'
            2' 3  3'
            3  3'
            ---------------
            1  2  2'
            2' 3  3'
            3  3'
        """        
        Type = self._type
        KingT = self.to_type('King')
        
        if i == 0:
            T0 = SemistandardTableau(KingT._rows)
            T1 = T0.bender_knuth_involution(1)
            Result = SymplecticTableau(self._n, rows = T1, type = 'King')
            return Result.to_type(Type)
        
        T0 = SemistandardTableau(KingT._rows)
        T1 = T0.bender_knuth_involution(2*i)
        T2 = T1.bender_knuth_involution(2*i-1)
        T3 = T2.bender_knuth_involution(2*i+1)
        T4 = T3.bender_knuth_involution(2*i)
        T5 = T4
        if len(T5) > i:
            while len(T5) < i or T5[i][0] == 2*i:
                D = {i:list(T5[i]) for i in range(len(T5))}
                row1 = list(T5[i-1][1:]) + [2*i+1]; row1.sort()
                row2 = list(T5[i][1:]) + [2*i+2]; row2.sort()
                D[i-1] = row1; D[i] = row2
                T5 = SemistandardTableau([D[i] for i in range(len(T5))])
        if display:
            print('Now applying Bender--Knuth to the tableau')
            tabs = [T0, T1, T2, T3, T4, T5]
            for i in range(len(tabs)):
                print('T%s :'%str(i))
                #print(GelfandTsetlinPattern(tabs[i]).pp())
                print(repr(SymplecticTableau(self._n, rows = tabs[i])))
                print('-'*15)
        Result = SymplecticTableau(self._n, rows = T5, type = 'King')
        return Result.to_type(Type)
    def to_GTpattern(self):
        r"""
        Returns the symplectic Gelfand--Tsetlin pattern associated with the tableau.
        
        EXAMPLES::
            sage: T = SymplecticTableau(3, rows = [[1,3,4],[5,5,6],[6,6]]); T
            1  2  2'
            3  3  3'
            3' 3'
            sage: T.to_GTpattern()
            3   3   2
            3   2   0
              3   0
              2   0
                1
                1
        """
        KingT = self.King()
        return KingTableau_to_SymplecticPattern(KingT)
    def King(self):
        r"""
        Returns the King tableau corresponding to the tableau.
        
        If ``self`` is a King tableau, this returns ``self``.
        
        If ``self`` is a DeConcini tableau, it applies Sheats' bijection.
        See :meth:`.SymplecticTableaux.Sheats`.
        
        If ``self`` is a Krattenthaler tableau, we first convert to a
        DeConcini tableau.
        
        If ``self`` is a Kashiwara tableau, we first convert to a
        Krattenthaler tableau.
        """
        if self._type == 'King':
            return self
        elif self._type == 'DeConcini':
            return DeConcini_to_King(self)
        elif self._type == 'Krattenthaler' or self._type == "split":
            return Krattenthaler_to_King(self)
        elif self._type == 'Kashiwara':
            return Kashiwara_to_King(self)
        raise ValueError("Cannot transfrom from custom to King tableau")
    def DeConcini(self):
        r"""
        Returns the De Concini tableau corresponding to the tableau.
        
        If ``self`` is a King tableau, it applies Sheats' bijection.
        See :meth:`.SymplecticTableaux.Sheats`.
        
        If ``self`` is a DeConcini tableau, this returns ``self``.
        
        If ``self`` is a Krattenthaler tableau, it takes the inverse
        of the split map. See :meth:`.SymplecticTableaux.splitInverse`.
        
        If ``self`` is a Kashiwara tableau, it first converts to a
        Krattenthaler tableau.
        """
        if self._type == 'King':
            return King_to_DeConcini(self)
        elif self._type == 'DeConcini':
            return self
        elif self._type == 'Krattenthaler' or self._type == "split":
            return Krattenthaler_to_DeConcini(self)
        elif self._type == 'Kashiwara':
            return Kashiwara_to_DeConcini(self)
        raise ValueError("Cannot transfrom from custom to De Concini tableau")
    def Krattenthaler(self):
        r"""
        Returns the Krattenthaler tableau corresponding to the tableau.
        
        If ``self`` is a King tableau, it first converts to a
        De Concini tableau.
        
        If ``self`` is a DeConcini tableau, it returns the split version.
        See :meth:`.SymplecticTableaux.splitVersion`.
        
        If ``self`` is a Krattenthaler tableau, this returns ``self``.
        
        If ``self`` is a Kashiwara tableau, it returns the cosplit version.
        See :meth:`.SymplecticTableaux.cosplitVersion`.
        """
        if self._type == 'King':
            return King_to_Krattenthaler(self)
        elif self._type == 'DeConcini':
            return DeConcini_to_Krattenthaler(self)
        elif self._type == 'Krattenthaler' or self._type == "split":
            return self
        elif self._type == 'Kashiwara':
            return Kashiwara_to_Krattenthaler(self)
        raise ValueError("Cannot transfrom from custom to Krattenthaler tableau")
    def Kashiwara(self):
        r"""
        Returns the Kashiwara tableau corresponding to the tableau.
        
        If ``self`` is a King tableau, it first converts to a
        De Concini tableau.
        
        If ``self`` is a DeConcini tableau, it first converts to a
        Krattenthaler tableau.
        
        If ``self`` is a Krattenthaler tableau, it takes the inverse
        of the cosplit map. See :meth:`.SymplecticTableaux.cosplitInverse`.
        
        If ``self`` is a Kashiwara tableau, this returns ``self``.
        """
        if self._type == 'King':
            return King_to_Kashiwara(self)
        elif self._type == 'DeConcini':
            return DeConcini_to_Kashiwara(self)
        elif self._type == 'Krattenthaler' or self._type == "split":
            return Krattenthaler_to_Kashiwara(self)
        elif self._type == 'Kashiwara':
            return self
        raise ValueError("Cannot transfrom from custom to Kashiwara tableau")
    def to_type(self, type):
        r"""
        Converts to the desired type.
        
        INPUT:
            - type: Either 'King', 'DeConcini', 'Krattenthaler', or 'Kashiwara'.
            
        EXAMPLES::
            sage: T = SymplecticTableau(3, rows = [[1,3],[5,5],[6,6]]); T
            1  2
            3  3
            3' 3'
            sage: T.to_type('DeConcini')
            3' 2'
            2' 1'
            2  3
            sage: T.to_type('Krattenthaler')
            1 1 2 2
            2 3 3 3
            4 5 6 6
            sage: T.to_type('Kashiwara')
            1  2
            3  3
            3' 1'
        """
        if type == 'King':
            return self.King()
        elif type == 'DeConcini':
            return self.DeConcini()
        elif type == 'Krattenthaler' or type == "split":
            return self.Krattenthaler()
        elif type == 'Kashiwara':
            return self.Kashiwara()
        raise ValueError("Can only convert to types 'King', 'DeConcini', 'Krattenthaler', or 'Kashiwara'")
    def f(self, i):
        r"""
        Applies the crystal operator f(i).
        
        Crystal operators are already implemented in Sage for Kashiwara's tableaux.
        See :meth:`crystals.tensor_product`. This functions converts the tableau
        to a Kashiwara tableau, performs the crystal operator, and converts back to
        the original type.
        
        """
        if self._type == 'King' or self._type == 'DeConcini':
            I = self._n+1-i
        elif self._type == 'Kashiwara' or self._type == 'Krattenthaler' or self._type == "split":
            I = i
        else:
            raise ValueError("Cannot perform crystal operators on custom tableaux")
            
        tabKashiwara = self.Kashiwara()
        tabCrystalElement = KashiwaraTableau_to_CrystalElement(tabKashiwara)
        tabCrystalElement = tabCrystalElement.f(I)
        if tabCrystalElement == None:
            return SymplecticTableau(self._n, rows = [], type = self._type)
        tabKashiwara = CrystalElement_to_KashiwaraTableau(tabCrystalElement, self._n)
        return tabKashiwara.to_type(self._type)
    def e(self, i):
        r"""
        Applies the crystal operator e(i).
        
        Crystal operators are already implemented in Sage for Kashiwara's tableaux.
        See :meth:`crystals.tensor_product`. This functions converts the tableau
        to a Kashiwara tableau, performs the crystal operator, and converts back to
        the original type.
        """
        if self._type == 'King' or self._type == 'DeConcini':
            I = self._n+1-i
        elif self._type == 'Kashiwara' or self._type == 'Krattenthaler' or self._type == "split":
            I = i
        else:
            raise ValueError("Cannot perform crystal operators on custom tableaux")
            
        tabKashiwara = self.Kashiwara()
        tabCrystalElement = KashiwaraTableau_to_CrystalElement(tabKashiwara)
        tabCrystalElement = tabCrystalElement.e(I)
        if tabCrystalElement == None:
            return SymplecticTableau(self._n, rows = [], type = self._type)
        tabKashiwara = CrystalElement_to_KashiwaraTableau(tabCrystalElement, self._n)
        return tabKashiwara.to_type(self._type)
        
class SymplecticTableauIterator:
    def __init__(self, tab):
        self._rows = tab._rows
        self._len = tab._len
        self._current_index = 0
    def __iter__(self):
        return self   
    def __next__(self):
        if self._current_index < self._len:
            self._current_index += 1
            return self._rows[self._current_index-1]        
        raise StopIteration
        
def transpose(tab):
    if tab == []:
        return []
    return [[row[i] for row in tab if len(row)>i] for i in range(len(tab[0]))]

def _circles(col, n):
    r'''
    Takes a column tableau (as a single array) in the `1<2<\cdots<2n`,
    relabels it to the alphabet `n'<\cdots<2'<1'<1<2<\cdots<n`, and returns
    its circle diagram as an array of two arrays.
    '''
    return [[n+1-i for i in col if i<=n]] + [[i-n for i in col if i>n]]
    
def splitVersion(tab):
    r'''
    Takes a tableau in the `1<2<\cdots<2n`, relabels it to
    `n'<\cdots<2'<1'<1<2<\cdots<n`, and returns its split version.
    '''
    new = []
    for col in tab._cols:
        new = new + _splitCol(col, tab._n)
    return SymplecticTableau(tab._n, cols = new, type = 'Krattenthaler')
    
def _splitCol(col, n):
    assert is_admissible(col, n), "this column is not admissible: %s"%(col)
    circ = _circles(col, n)
    A = circ[0]
    D = circ[1]
    H = [i for i in [1..n] if i not in A and i not in D]
    I = [i for i in [1..n] if i     in A and i     in D]
    J = []; H_ = H
    for i in [1..len(I)]:
        j = max(h for h in H_ if h < I[i-1])
        H_.remove(j); J.append(j)
    B = [i for i in [1..n] if (i in A and i not in I) or i in J]; B.reverse()
    C = [i for i in [1..n] if (i in D and i not in I) or i in J]
    return [[n+1-i for i in A] + [n+i for i in C], [n+1-i for i in B] + [n+i for i in D]]
    
def cosplitVersion(tab):
    r'''
    Takes a tableau in the `1<2<\cdots<2n`, relabels it to
    `n'<\cdots<2'<1'<1<2\cdots<n`, and returns its cosplit version.
    '''
    new = []
    for col in tab._cols:
        new = new + _cosplitCol(col, tab._n)
    return SymplecticTableau(tab._n, cols = new, type = 'Krattenthaler')
    
def _cosplitCol(col, n):
    assert is_coadmissible(col, n), "this column is not coadmissible: %s"%(col)
    circ = _circles(col, n)
    B = circ[0]
    C = circ[1]
    H = [i for i in [1..n] if i not in B and i not in C]
    I = [i for i in [1..n] if i     in B and i     in C]
    J = []; H_ = H
    for i in [1..len(I)]:
        j = min(h for h in H_ if h > I[i-1])
        H_.remove(j); J.append(j)
    A = [i for i in [1..n] if (i in B and i not in I) or i in J]; A.reverse()
    D = [i for i in [1..n] if (i in C and i not in I) or i in J]
    return [[n+1-i for i in A] + [n+i for i in C], [n+1-i for i in B] + [n+i for i in D]]
    
def splitInverse(tab):
    r'''
    Takes the cosplit version of a tableau (in the alphabet
    `1<2<\cdots<2n`) and gives the original tableau.    
    '''
    new = _splitInverseCols(tab._cols, tab._n)
    return SymplecticTableau(tab._n, cols = new, type = 'DeConcini')

def _splitInverseCols(cols, n):
    new = []
    for i in range(len(cols)/2):
        new.append([j for j in cols[2*i] if j <= n] + [j for j in cols[2*i+1] if j>n]) 
    return new
    
def cosplitInverse(tab):
    r'''
    Takes the cosplit version of a tableau (in the alphabet
    `1<2<\cdots<2n`) and gives the original tableau.    
    '''    
    new = _cosplitInverseCols(tab._cols, tab._n)
    return SymplecticTableau(tab._n, cols = new, type = 'Kashiwara')
    
def _cosplitInverseCols(cols, n):
    new = []
    for i in range(len(cols)/2):
        new.append([j for j in cols[2*i+1] if j <= n] + [j for j in cols[2*i] if j>n]) 
    return new
        
def is_admissible(c, n):
    r'''
    Takes a semistard column in `1<2<\cdots<2n`, relabels it to
    `n'<\cdots<2'<1'<1<2<\cdots<n` and checks admissibility.
    '''
    for i in [1..n]:
        iPrime = n+1-i
        iNormal = n+i
        entries = [e for e in c if e >= iPrime and e <= iNormal]
        if len(entries) > i:
            return False
    return True  
    
def is_coadmissible(c, n):
    r'''
    Takes a semistard column in `1<2<\cdots<2n`, relabels it to
    `n'<\cdots<2'<1'<1<2<\cdots<n` and checks coadmissibility.
    '''
    for i in [1..n]:
        iPrime = n+1-i
        iNormal = n+i
        entries = [e for e in c if e <= iPrime or e >= iNormal]
        if len(entries) > n+1-i:
            return False
    return True

def King_to_DeConcini(tab):
    return SheatsInverse(tab)
def DeConcini_to_King(tab):
    return Sheats(tab)
def King_to_Kashiwara(tab):
    return DeConcini_to_Kashiwara(King_to_DeConcini(tab))
def Kashiwara_to_King(tab):
    return DeConcini_to_King(Kashiwara_to_DeConcini(tab))
def DeConcini_to_Kashiwara(tab):
    return cosplitInverse(splitVersion(tab))
def Kashiwara_to_DeConcini(tab):
    return splitInverse(cosplitVersion(tab))
def King_to_Krattenthaler(tab):
    return splitVersion(King_to_DeConcini(tab))
def Krattenthaler_to_King(tab):
    return DeConcini_to_King(splitInverse(tab))
def DeConcini_to_Krattenthaler(tab):
    return splitVersion(tab)
def Kashiwara_to_Krattenthaler(tab):
    return cosplitVersion(tab)
def KashiwaraTableau_to_CrystalElement(tab):
    assert tab._type == 'Kashiwara', "This function expects a Kashiwara tableau"
    import sage.combinat.crystals.tensor_product
    Tab = crystals.Tableaux(['C',tab._n], shape = tab._shape)
    rows = [[i if i<=tab._n else i-2*tab._n-1 for i in row] for row in tab._rows]
    return Tab(rows = rows)
def CrystalElement_to_KashiwaraTableau(tab, n):
    rows = [[i if i>=0 else i+2*n+1 for i in row] for row in tab.to_tableau()]
    return SymplecticTableau(n, rows = rows, type = 'Kashiwara')
    
#############################################################################################################
def Sheats(tab):
    r"""
    Given a tableau in the alphabet `1<2<\cdots<2n`, it is
    interpreted as a De Concini tableau in the alphabet
    `n'<\cdots<2'<1'<1<2<\cdots<n`. Then, the Sheats algorithm
    is performed on the tableau, resulting on a King tableau
    in the alphabet `1<2<\cdots<2n` (interpreted as in the 
    alphabet `1'<1<\cdots<n'<n`).
    
    EXAMPLES::
        sage: tab = SymplecticTableau(3, rows = [[1, 2], [2, 3], [5, 6]], type = 'DeConcini'); tab
        3' 2'
        2' 1'
        2  3
        sage: Sheats(tab)
        1  2
        3  3
        3' 3'
        sage: tab = SymplecticTableau(3, rows = [[1,3,3], [2,5], [4,6]], type = 'DeConcini'); tab
        3' 1' 1'
        2' 2
        1  3
        sage: Sheats(tab)
        1  2  2
        2' 2'
        3  3'
    """
    # we will be perfoming transformations on a subtableau and then
    # adding some stuff to each column.
    lam1 = tab._shape[0]
    add = [[] for i in range(lam1)]
    n = tab._n
    cols = tab._cols

    # n gives the size of the alphabet
    for m in [1..n]:
        # m is the entry to apply jdtq on
        # m is taken from the alphabet 1<2<\cdots<2n
        mPrime = m
        mNormal = 2*n+1-m
        
        # in the new order, the entries that we moved have another relative position   
        newMPrime = 2*(n-m)+1
        newMNormal = 2*(mNormal-n)
        
        # k is the number of times m' appears
        k = len([col for col in cols if len(col)!=0 and col[0] == mPrime])
        
        # the entries equal to m will stay fixed.
        # but in the new tableau, they are in a different
        # position of the alphabet, so translation is needed.
        colsFix = [[newMNormal for i in col if i == mNormal] for col in cols]
        colsFix = colsFix + [[] for i in range(lam1 - len(colsFix))]
        # everything else that is not m or m' will be movable
        colsMove = [[i for i in col if i != mNormal and i != mPrime] for col in cols]
        
        # we will now transform colsMove and later add some extra entries
        for i in range(lam1):
            add[i] = colsFix[i] + add[i]
        # this transformation expects a split tableau
        Cols = sum([_splitCol(col, n) for col in colsMove], [])
        
        while k != 0:
            # we apply _sjdt to the only inner corner
            (Cols, (p,q), n) = _sjdt(Cols, (1,k), n)
            
            # we keep track of the entry we need to add
            # later to the tableau
            cols = _splitInverseCols(Cols, n)
            cols = cols + [[] for i in range(lam1-len(cols))]
            add[q-1] = [newMPrime] + add[q-1]
            for i in range(lam1):
                add[i] = [newMNormal]*cols[i].count(mNormal) + add[i]
            
            # we prepare for the next iteration
            k = k + len([col for col in cols if mPrime in col]) - 1
            colsMove = [[i for i in col if i != mNormal and i != mPrime] for col in cols]
            Cols = sum([_splitCol(col, n) for col in colsMove], [])
        colsMoved = _splitInverseCols(Cols, n)
        colsMoved = colsMoved + [[] for i in range(lam1- len(colsMoved))]
        # colsMoved = [[i-1 for i in col] for col in colsMoved]
         
        # we prepare the body for the next iteration
        cols = colsMoved
    # now comes the time to add the entries that we were keeping aside
    newCols = []
    for i in range(lam1):
        newCols.append(colsMoved[i] + add[i])
    
    return SymplecticTableau(tab._n, cols = newCols, type = 'King')

def _sjdt(cols, puncture, n):
    r"""
    Takes a set of columns that form a punctured split skew tableau
    if the puncture is set at the input location. We let 0s mark the
    empty skew boxes.
    
    The puncture is giving coordinates in terms of the non-split
    version of the tableau, starting at (1,1).
    """
    (p,q) = puncture # starts at (1,1)
    
    # Step 1:
    # if the puncture is an outer corner, stop.
    if len(cols) <= 2*q or len(cols[2*q]) < p:
        if len(cols[2*q-1]) < p:
            return (cols, (p,q), n)
    # Step 2:
    # if the column to the right is too short,
    # just move the puncture down.
        return _sjdt(cols, (p+1, q), n)
    # Step 3:
    # if the current column is too short,
    # move the puncture to the right.
    if len(cols[2*q-1]) < p:
        return _sjdt_subroutine(cols, puncture, n)
    # Step 4:
    # when the algorithm reaches this step, it
    # has checked that there are two possible
    # locations to move the puncture to.
    # We now compare entries to see where to go.
    if cols[2*q-1][p-1] <= cols[2*q][p-1]:
        return _sjdt(cols, (p+1, q), n)
    else:
        return _sjdt_subroutine(cols, puncture, n)

def _sjdt_subroutine(cols, puncture, n):
    (p, q) = puncture
    a = cols[2*q][p-1]
    
    if a <= n: # if left is primed
        leftBuffer = cols[2*q-2].count(0)
        leftCurrent = [cols[2*q-2][leftBuffer:], cols[2*q-1][leftBuffer:]]
        left = _cosplitInverseCols(leftCurrent, n)[0]
        leftCircs = _circles(left, n)
        # leftCircs = [[i for i in leftCircs[0] if i>a] + [n+1-a] + [i for i in leftCircs[0] if i<a], leftCircs[1]]
        # leftArray = [n+1-i for i in leftCircs[0]] + [i+n for i in leftCircs[1]]
        leftArray = [n+1-i for i in leftCircs[0]] \
                    + [a] \
                    + [i+n for i in leftCircs[1]]
        leftArray.sort()           
        left = _cosplitCol(leftArray, n)
        left = [[0]*leftBuffer + list(left[0]), [0]*leftBuffer + list(left[1])]
        
        rightBuffer = cols[2*q].count(0)
        rightCurrent = [cols[2*q][rightBuffer:], cols[2*q+1][rightBuffer:]]
        right = _splitInverseCols(rightCurrent, n)[0]
        rightCircs = _circles(right, n)
        rightCircs[0].remove(n+1-a)
        rightArray = [n+1-i for i in rightCircs[0]] + [i+n for i in rightCircs[1]]
        right = _splitCol(rightArray, n)
        right = [[0]*rightBuffer + list(right[0]), [0]*rightBuffer + list(right[1])]
        
        return _sjdt(cols[:2*q-2] + left + right + cols[2*q+2:], (p, q+1), n)
    else: # if left is non-primed
        leftBuffer = cols[2*q-2].count(0)
        leftCurrent = [cols[2*q-2][leftBuffer:], cols[2*q-1][leftBuffer:]]
        left = _splitInverseCols(leftCurrent, n)[0]
        leftCircs = _circles(left, n)
        # leftCircs = [leftCircs[0], [i for i in leftCircs[1] if i<a] + [a-n] + [i for i in leftCircs[1] if i>a]]
        # leftArray = [n+1-i for i in leftCircs[0]] + [i+n for i in leftCircs[1]]
        leftArray = [n+1-i for i in leftCircs[0]] \
                    + [a] \
                    + [i+n for i in leftCircs[1]]
        leftArray.sort()
        left = _splitCol(leftArray, n)
        left = [[0]*leftBuffer + list(left[0]), [0]*leftBuffer + list(left[1])]
        
        rightBuffer = cols[2*q].count(0)
        rightCurrent = [cols[2*q][rightBuffer:], cols[2*q+1][rightBuffer:]]
        right = _cosplitInverseCols(rightCurrent, n)[0]
        rightCircs = _circles(right, n)
        rightCircs[1].remove(a-n)
        rightArray = [n+1-i for i in rightCircs[0]] + [i+n for i in rightCircs[1]]
        right = _cosplitCol(rightArray, n)
        right = [[0]*rightBuffer + list(right[0]), [0]*rightBuffer + list(right[1])]
        
        return _sjdt(cols[:2*q-2] + left + right + cols[2*q+2:], (p, q+1), n)

def SheatsInverse(tab):
    r"""
    Given a tableau in the alphabet `1<2<\cdots<2n`, it is
    interpreted as a King tableau in the alphabet
    `1'<1<2'<2<\cdots<n'<n`. Then, the Sheats inverse algorithm
    is performed on the tableau, resulting on a De Concini tableau
    in the alphabet `1<2<\cdots<2n` (interpreted as in the 
    alphabet `n'<\cdots<1'<1<\cdots<n`).
    
    EXAMPLES::
        sage: tab = SymplecticTableau(3, rows = [[1, 3], [5, 5], [6, 6]]); tab
        1  2
        3  3
        3' 3'
        sage: SheatsInverse(tab)
        3' 2'
        2' 1'
        2  3
        sage: tab = SymplecticTableau(3, rows = [[1, 3, 3], [4, 4], [5, 6]]); tab
        1  2  2
        2' 2'
        3  3'
        sage: SheatsInverse(tab)
        3' 1' 1'
        2' 2
        1  3
    """
    cols = tab._cols
    n = tab._n
    # alphabet = sum([[str(i) + "'", str(i) + " "] for i in [1..n]], [])
    for m in [2..n]:
        mPrime = 2*m-1
        mNormal = 2*m
        newMPrime = 1
        newMNormal = 2*m
        # newalphabet = [alphabet[mPrime-1]] + alphabet[:mPrime-1] + alphabet[mPrime:]
        
        colsMove = [[i+1 for i in col if i < mPrime] for col in cols]
        punctures = sum([[(i+1, j+1)
                    for i in range(len(cols[j])) 
                    if cols[j][i] == mPrime]
                    for j in range(len(cols))], [])
        colsFix = [[i for i in col if i > mPrime] for col in cols]
        
        while len(punctures) != 0:
            puncture = punctures[0]; punctures = punctures[1:]
            
            Cols = sum([_splitCol(col, m) for col in colsMove], [])
            
            (ColsSig, puncSig, _) = _sigmaAux(Cols, puncture, m)
            (ColsSig, puncSig, _) = _sjdt(ColsSig, puncSig, m)
            (Cols, (_, _), _) =_sigmaAux(ColsSig, puncSig, m)
            
            colsMove = _splitInverseCols(Cols, m)
                        
            k = len([col for col in colsMove if len(col)!=0 and col[0]==newMPrime])
            for col in colsMove:
                while 0 in col:
                    col.remove(0)
            if mNormal in colsFix[k]:
                colsMove[k] = colsMove[k] + [newMNormal]
                colsFix[k].remove(mNormal)
            colsMove[0] = [newMPrime] + colsMove[0]
            
        cols = [colsMove[i] + colsFix[i] for i in range(len(cols))]
        for col in cols:
            while 0 in col:
                col.remove(0)
        # alphabet = newalphabet
    return SymplecticTableau(n, cols = cols, type = 'DeConcini')

def _sigmaAux(cols, puncture, n):
    r"""
    This function does several things:
      - Rotates a skew split tableau 180 degrees.
      - Primes every non-primed entry and unprimes every primed entry.
      - Gives the corresponding puncture location after rotation.
    """
    (p, q) = puncture
    
    lam1 = len(cols)
    length = len(cols[0])
    
    Cols = [[2*n+1-col[i] if (i < len(col) and col[i] != 0) else 0 for i in range(length)] for col in cols]
    if 0 in Cols[2*q-1]:
        Cols[2*q-1].remove(0)
    if 0 in Cols[2*q-2]:
        Cols[2*q-2].remove(0)
    
    new = []
    for i in range(lam1)[::-1]:
        col = deepcopy(Cols[i])
        col.reverse()
        new.append(col)
        
    newPuncture = (length - p + 1, int(lam1/2 - q + 1))
    if newPuncture[0] == 0:
        newPuncture = (newPuncture[0]+1, newPuncture[1])
        new = [[0] + col for col in new]
    return (new, newPuncture, n)

#############################################################################################################
    
def check_all_BKinvolutions(lam, i, max_entry = None):
    if max_entry == None:
        max_entry = 2*len(lam)
    for t in SemistandardTableaux(lam, max_entry = max_entry):
        t = SymplecticTableau(max_entry, rows = t)
        if t.is_well_defined():
            if t != t.bender_knuth_involution(i).bender_knuth_involution(i):
                print(t)
                break
    print('Test passed!')
    
class SymplecticPattern(SageObject):
    def __init__(self, input):
        self._list = input
        self._len = len(input)
    def list(self):
        return self._list
    def len(self):
        return self._len
    def _repr_(self):
        return '\n'.join(''.join(['  ']*floor(i/2) + ['   '.join(map(str,self._list[i]))]) for i in range(self._len))
    def to_tableau(self):
        return SymplecticPattern_to_KingTableau(self)

def KingTableau_to_SymplecticPattern(T):
    L = T.list(); N = 2*T._n
    G = GelfandTsetlinPattern(Tableau(L))
    return SymplecticPattern([G[i][:int(N/2 - floor(i/2))] for i in range(N)])
    
def SymplecticPattern_to_KingTableau(G):
    """
    EXAMPLES::
        sage: G = SymplecticPattern([[3,2,0],[3,0,0],[3,0],[3,0],[3],[1]]); G
        3   2   0
        3   0   0
          3   0
          3   0
            3
            1
        sage: G.to_tableau()
        1  1' 1'
        3' 3'
    """
    L = G.list(); n = G.len()
    D = {cell : n for cell in Partition(L[0]).cells()}
    for (Part, i) in zip(L, range(n,0,-1)):
        for cell in Partition(Part).cells():
            D[cell] = i
    return SymplecticTableau(n, rows = _dict_to_tableau(D, L[0]))
    
def _dict_to_tableau(D, shape):
    rows = []
    for i in range(len(shape)):
        row = []
        for j in range(shape[i]):
            row.append(D[(i,j)])
        if row != []:
            rows.append(row)
    return Tableau(rows)
    
    
class abstractPattern(SageObject):
    def __init__(self, dict = None):
        n = 3; self._n = n
        if dict != None:
            self._dict = dict
        else:
            self._dict = {
                j : {
                    i : var('x%s%s'%(i,j))
                    for i in [1..ceil(j/2)]
                }
                for j in [1..2*n][::-1]
            }
            self._dict[4][3] = 0
    def dict(self):
        return self._dict
    def list(self):
        return [[self._dict[rowKey][colKey] for colKey in self._dict[rowKey]] for rowKey in self._dict]
    def size(self):
        return self._n
    def _repr_(self):
        return self.list()
    def toggle(i, j):
        D = self._dict
        if j == 2*i:
            D[j][j] = 0
        D[i][j] = Min(D[j+1][i], D[j-1][i-1]) + Max(D[j+1][i+1], D[j-1][i]) - D[j][i]
        return abstractPattern(dict = D)
