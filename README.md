# Projection via \(Q Q^T\)

This document describes the matrix transformation 
\[
P = Q Q^T,\quad P Q = (Q Q^T) Q
\]
where \(Q\) is an \(n \times m\) matrix with orthonormal columns (\(Q^T Q = I_m\)). Its purpose is to formalize and document the key properties and operations of the orthogonal projection defined by \(Q\).

---

## 🧮 Definitions

- **Matrix**: \(Q \in \mathbb{R}^{n \times m}\)
- **Orthonormality**: \(Q^T Q = I_m\)
- **Projection Operator**:  
  \[
    P = Q Q^T
  \]

---

## Properties of \(P = Q Q^T\)

1. **Symmetry**  
   \[
     P^T = (Q Q^T)^T = Q Q^T = P
   \]

2. **Idempotence**  
   \[
     P^2 = Q Q^T Q Q^T = Q (Q^T Q) Q^T = Q I_m Q^T = P
   \]

3. **Orthogonal Projection**  
   \(P\) projects any vector \(x \in \mathbb{R}^n\) onto the column space of \(Q\).  
   - If \(x\in \mathrm{col}(Q)\), then \(P x = x\).  
   - If \(x\) is orthogonal to \(\mathrm{col}(Q)\), then \(P x = 0\).  
   :contentReference[oaicite:1]{index=1}

4. **Action on \(Q\) Itself**  
   \[
     (Q Q^T) Q = Q (Q^T Q) = Q I_m = Q
   \]
   Thus, columns of \(Q\) are fixed by this projection.

5. **Eigenstructure**  
   - Eigenvalue \(1\): multiplicity \(m\), eigenvectors = columns of \(Q\).  
   - Eigenvalue \(0\): multiplicity \(n - m\), eigenvectors orthogonal to \(\mathrm{col}(Q)\).

---

## Geometric Interpretation

- **Shape preservation**: \(P\) preserves vectors already in the column space of \(Q\), and annihilates those in its orthogonal complement.  
- **Subspace invariance**: \(\mathrm{col}(Q)\) is invariant under \(P\); its orthogonal complement is the null space of \(P\).  
  :contentReference[oaicite:2]{index=2}

---

## Example

Let 
```text
Q = [u₁ | ··· | uₘ], where {uᵢ} is an orthonormal basis of a subspace U ⊆ ℝⁿ.
