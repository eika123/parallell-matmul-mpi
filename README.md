# Distributed matrix multiplication and eigenvalue/eigenvector calculations
Simple distributed matrix multiplications using MPI, with power iteration to compute the principal eigenvector.
The principal eigenvector corresponds to the largest eigenvalue.

In it's current state just makes random matrices and computes eigenvectors.

At least back in the day, power iteration was at the foundation of googles [PageRank algorithm](https://en.wikipedia.org/wiki/PageRank).
