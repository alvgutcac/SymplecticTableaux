# SymplecticTableaux
A SageMath library implementing symplectic tableaux and symplectic patterns.

### Instalation
To use this locally, save the file in your computer (e.g. in Desktop). Then run a Sage console in this path directory and load the file.
```
sage: cd Desktop
/home/sage/Desktop/
sage: load("SymplecticTableaux.sage")
```

### Documentation
Documentation will be generated at a later stage. For now, some functions have some preliminary documentation as docstrings. This can be accessed also from the console.
```
sage: SymplecticTableau.f?
Signature:      SymplecticTableau.f(self, i)
Docstring:        Applies the crystal operator f(i).
        
                  Crystal operators are already implemented in Sage for Kashiwara's tableaux.
                  See :meth:`crystals.tensor_product`. This functions converts the tableau
                  to a Kashiwara tableau, performs the crystal operator, and converts back to
                  the original type.
Init docstring: Initialize self.  See help(type(self)) for accurate signature.
File:           Dynamically generated function. No source code available.
Type:           function
```
