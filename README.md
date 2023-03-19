# Communication-efficient-secret-sharing

Skoltech Intro to Blockchain final project

Topic: Communication efficient secret sharing 
  (https://arxiv.org/pdf/1505.07515.pdf)

Team Members:\
\
Andreev Nikita\
Begichev Ilya\
Galichin Andrey

About:

In this project, we consider the implementation of Shamir secret scheme, and it's improvement - https://arxiv.org/pdf/1505.07515.pdf. First one could be found in shamir.py script, latter one - in efficient.py. We provide a simple class-based functional to operate with both of these approaches.\
An example on how to use our code can be seen in example.ipynb notebook.

Requirements:

As our main building block we use $\textbf{galois}$ python package.\
galois library is a Python 3 package that extends NumPy arrays to operate over finite fields.\
This is particular use in our case, so we use it to efficiently operate with polynomials.\

Download:
```console
  pip install galois
```
