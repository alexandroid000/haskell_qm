---
title: Quantum Vector
author: Jan Skibinski
...

Quantum vector
==============

**Jan Skibinski, [Numeric Quest Inc.](http://www.numeric-quest.com/news/),
Huntsville, Ontario, Canada

Literate Haskell module *QuantumVector.lhs***

Initialized: 2000-05-31, last modified: 2000-06-10

Introduction
------------

This is our attempt to model the abstract Dirac's formalism of Quantum Mechanics
in Haskell. Although we have been developing quantum mechanical applications and
examples for some time [2], the machinery used there is tightly coupled to a
concrete representation of states and observables by complex vectors and
matrices. implemented mainly as Haskell lazy lists.

However, the Dirac's formalism in Hilbert space is much more abstract than that,
and many problems of Quantum Mechanics can be solved without referring to any
particular matrix representation, but using certain generic properties of
operators, such as their commutative relations instead. Haskell seems to be well
suited for such abstract tasks, even in its current form that does not support
any of the abstract notions of computer algebra as yet. This has been already
recognized by Jerzy Karczmarczuk [1], where he proposes a very interesting
representation of Hilbert space and illustrates it by several powerful examples.
But the task is not trivial and far from being complete. Quantum Mechanics
presents many challenges to any formalism and only by careful examination of
many of its facets and alternative approaches, a consistent model of Dirac's
formalism can be developed for Haskell. Hoping to help with solving this
problem, we present here a computing abstract, which is quite different from
that of [1].

We recognize a quantum state as an abstract vector | x \>, which can be
represented in one of many possible bases -- similar to many alternative
representations of a 3D vector in rotated systems of coordinates. A choice of a
particular basis is controlled by a generic type variable, which can be any
Haskell object -- providing that it supports a notion of equality and ordering.
A state which is composed of many quantum subsystems, not necessarily of the
same type, can be represented in a vector space considered to be a tensor
product of the subspaces.

With this abstract notion we proceed with Haskell definition of two vector
spaces: Ket and its dual Bra. We demonstrate that both are properly defined
according to the abstract mathematical definition of vector spaces. We then
introduce inner product and show that our Bra and Ket can be indeed considered
the vector spaces with inner product.  Multitude of examples is attached in the
description. To verify the abstract machinery developed here we also provide the
basic library module
[Momenta](http://www.numeric-quest.com/haskell/Momenta.html) -- a non-trivial
example designed to compute Clebsch-Gordan coefficients of a transformation from
one basis of angular momenta to another.

Section 6 is a rehash of known definitions of linear operators with the emphasis
on both Dirac and Haskell notations and on Haskell examples. The formalism
developed here centers around two operations: a scalar product of two vectors,
**x \<\> y**, and a closure operation, **a \>\< x**, which can be considered an
application of a quantum operator **a** to a vector **x**. At this stage our
formalism applies only to discrete cases, but we hope to generalize it on true
Hilbert space as well.

Contents
--------

-   1. Infix operators
-   2. Vector space
-   3. Ket vector space
-   4. Bra vector space
-   5. Bra and Ket spaces as inner product spaces
-   6. Linear operators
    -   6.1. Operator notation
    -   6.2. Renaming the representation
    -   6.3. Closure formula, or identity operator
    -   6.4. Changing the representation
    -   6.5. Implementation of the operator equation A | x \> = | y \>
    -   6.6. Inverse operator
    -   6.7. Matrix representation of an operator
    -   6.8. Adjoint operator
    -   6.9. Unitary operator
    -   6.10. Hermitian operator
-   7. Showing kets and bras
-   8. Data Tuple for tensor products
-   9. References
-   10. Copyright and license


Infix operators
===============

Haskell requires that fixities of infix operators are defined at the top of the
module. So here they are. They are to be explained later.

> module QuantumVector where
> import Data.Complex -- our Scalar is Complex Double
> import Data.List (nub)

> infixl 7 *>  -- tensor product of two kets
> infixl 7 <*  -- tensor product of two bras

> infix 6 |>   -- scalar-ket multiplication
> infix 6 <|   -- scalar-bra multiplication

> infixl 5 +>  -- sum of two kets
> infixl 5 <+  -- sum of two bras

> infix 4 <>  -- inner product
> infix 5 ><  -- closure

> infix 5 >-<
> infix 5 >+<
> infix 5 >*<

infix 5 >*
infix 4 *<


Vector space
============

Definition. A set V of elements x ,y ,z ,...is called a vector (or linear) space
over a complex field C if

-   vector addition + is defined in V such that V is an abelian group under
    addition, with identity element 0

1: x + y       = y + x
2: x + (y + z) = (x + y) + z
3: 0 + x       = x + 0

-   the set is close with respect to scalar multiplication and vector addition

4: a (x + y)   = a x + a y
5: (a + b) x   = a x + b x
6: a (b x)     = (a b) x
7: 1 x         = x
8: 0 x         = 0
    where
        a, b, c are complex scalars

Definition. The maximum number of linearly independent vectors in V or, what is
the same thing, the minimum number of linearly independent vectors required to
span V is the dimension r of vector space V.

Definition. A set of r linearly independent vectors is called a basis of the
space. Each vector of the space is then a unique linear combination of the
vectors of this basis.

Based on the above definitions we will define two vector spaces: ket space and
its dual -- bra space, which, in addition to the above properties, will also
support several common operations -- grouped below in the class DiracVector.

> class DiracVector a where
>     zero       :: a
>     add        :: a -> a -> a
>     scale      :: Scalar -> a -> a
>     reduce     :: a -> a
>     basis      :: a -> [a]
>     components :: a -> [Scalar]
>     compose    :: [Scalar] -> [a] -> a
>     dimension  :: a -> Int
>     norm       :: a -> Double
>     normalize  :: a -> a

>     dimension x   = length (basis x)

>     normalize x
>         | normx == 0 = x
>         | otherwise  = compose cs (basis x)
>          where
>             cs     = [a*v :+ b*v |a :+ b <- components x]
>             v      = 1 / normx
>             normx  = norm x


Ket vector space
================

We submit that the following datatype and accompanying operations define a
complex vector space, which we will call the ket vector space.

> type Scalar = Complex Double

> data Ket a  =
>            KetZero                     -- zero ket vector
>          | Ket a                       -- base ket vector
>          | Scalar  :|> Ket a           -- scaling ket vectors
>          | Ket a   :+> Ket a           -- spanning ket space

A tensor product of two ket spaces is also a ket space.

> (*>) :: (Ord a, Ord b) => Ket a -> Ket b -> Ket (Tuple a b)
> Ket a   *> Ket b    = Ket (a :* b)
> _       *> KetZero  = KetZero
> KetZero *> _        = KetZero
> x       *> y        = foldl1 (:+>) [((Bra a <> x) * (Bra b <> y)) :|> Ket (a :* b)
>                                   | Ket a <- basis x, Ket b <- basis y]

> (|>) :: Ord a => Scalar -> Ket a -> Ket a
>     --
>     -- Multiplication of ket by scalar
>     --
> s |> (x :+> y)  = (s |> x) +> (s |> y)
> _ |> KetZero    = KetZero
> 0 |> _          = KetZero
> s |> (s2 :|> x) = (s * s2) |> x
> s |> x          = s :|> x


> (+>) :: Ord a => Ket a  -> Ket a  -> Ket a
>     --
>     -- Addition of two kets
>     --
> x +> KetZero = x
> KetZero +> x = x
> x +> y       = reduce (x :+> y)


> instance (Eq a, Ord a) => Eq (Ket a) where
>     --
>     -- Two ket vectors are equal if they have identical
>     -- components
>     --
>     x == y = and [c k x == c k y  | k <- basis x]
>         where
>             c k z = (toBra k) <> z

The data Ket is parametrized by type variable "a", which can be anything that
can be compared for equality and ordered: integer, tuple, list of integers, etc.
For example, the data constructor `Ket (3::Int)` creates a base vector `|3>`,
annotated by Int.  Similarly, `Ket (2::Int,1::Int)`, creates a base vector
`|(2,1)>` annotated by a tuple of Ints. Those two vectors belong to two
different bases.

The eight examples below illustrate the eight defining equations of the vector
space, given in section 1. All of them evaluate to True.

ghci> Ket 2 +> Ket 3            == Ket 3 +> Ket 2
ghci> Ket 1 +> (Ket 2 +> Ket 3) == (Ket 1 +> Ket 2) +> Ket 3
ghci> Ket 1 +> KetZero          == KetZero +> Ket 1
ghci> 5 |> (Ket 2 +> Ket 3)     == 5 |> Ket 2 +> 5 |> Ket 3
ghci> (5 + 7) |> Ket 2          == 5 |> Ket 2 +> 7 |> Ket 2
ghci> 2 |> (4 |> Ket 2)         == 8 |> Ket 2
ghci> 1 |> Ket 2                == Ket 2
ghci> 0 |> Ket 2                == KetZero

The ket expressions can be pretty printed, as shown below.

ghci> Ket 2 +> Ket 3        ==> 1.0 |2> + 1.0 |3>
ghci> 5 |> (Ket 2 +> Ket 3) ==> 5.0 |2> + 5.0 |3>
ghci> 2 |> (4 |> Ket 2)     ==> 8.0 |2>

In order to support all those identities we also need several additional
functions for reducing the vector to its canonical form, for composing the ket
vector, and for extracting the ket basis and the ket components -- as shown
below.

> reduceKet :: Ord a => Ket a -> Ket a
> reduceKet x
>     --
>     -- Reduce vector `x' to its canonical form
>     --
>     = compose cs ks
>       where
>           ks = basis x
>           cs = [toBra k <> x | k <- ks]

> ketBasis :: Ord a => Ket a -> [Ket a]
>     --
>     -- Sorted list of unique base vectors of the ket vector
>     --
> ketBasis KetZero        = []
> ketBasis (Ket k)        = [Ket k]
> ketBasis (_ :|> x)      = [x]
> ketBasis (k1 :+> k2)    = nub (ketBasis k1 ++ ketBasis k2)

> toBra :: Ord a => Ket a -> Bra a
>     --
>     -- Convert from ket to bra vector
>     --
> toBra (Ket k)           = Bra k
> toBra (x :+> y)         = toBra x :<+ toBra y
> toBra (p :|> x)         = (conjugate p) :<| toBra x

> instance Ord a => DiracVector (Ket a)  where
>     zero          = KetZero
>     add           = (+>)
>     scale         = (|>)
>     reduce        = reduceKet
>     basis         = ketBasis
>     components x  = [toBra e <> x | e <- basis x]
>     compose xs ks = foldl1 (:+>) [fst z :|> snd z | z <- zip xs ks]

>     norm KetZero  = 0
>     norm x        = sqrt $ realPart (toBra x <> x)

But those auxilliary functions refer to vectors from the conjugated space bra,
which we shall now define below.


Bra vector space
================

Definition. Let V be the defining n-dimensional complex vector space.  Associate
with the defining n-dimensional complex vector space V a conjugate (or dual)
n-dimensional vector space obtained by complex conjugation of elements x in V.

We will call this space the bra space, and the corresponding vectors - the bra
vectors. Further, we submit that the following datatype and the corresponding
operations define bra space in Haskell.

> data Bra a =
>            BraZero                   -- zero bra vector
>          | Bra a                     -- base bra vector
>          | Scalar :<| Bra a          -- scaling bra vectors
>          | Bra a  :<+ Bra a          -- spanning bra space

A tensor product of two bra spaces is also a bra space.

> (<*) :: (Ord a, Ord b) => Bra a -> Bra b -> Bra (Tuple a b)
> Bra a   <* Bra b    = Bra (a :* b)
> _       <* BraZero  = BraZero
> BraZero <* _        = BraZero
> x       <* y        = foldl1 (:<+) [((x <> Ket a) * (y <> Ket b)) :<| Bra (a :* b)
>                                   | Bra a <- basis x, Bra b <- basis y]

> (<|) :: Ord a => Scalar -> Bra a -> Bra a
> s <| (x :<+ y)  = (s <| x) <+ (s <| y)
> _ <| BraZero    = BraZero
> 0 <| _          = BraZero
> s <| (s2 :<| x) = (s * s2) <| x
> s <| x          = s :<| x

> (<+) :: Ord a => Bra a -> Bra a -> Bra a
>     --
>     -- Sum of two bra vectors
>     --
> x <+ BraZero = x
> BraZero <+ x  = x
> x <+ y       = reduce (x :<+ y)

> instance (Eq a, Ord a) => Eq (Bra a) where
>     --
>     -- Two bra vectors are equal if they have
>     -- identical components
>     --
>     x == y = and [c b x == c b y  | b <- basis x]
>         where
>             c b z = z <> toKet b

Similarly to what we have done for ket vectors, we also define several
additional functions for reducing the bra vector to its canonical form, for
composing the bra vector, and for extracting the bra basis and the bra
components -- as shown below.

> reduceBra :: Ord a => Bra a -> Bra a
> reduceBra x
>     --
>     -- Reduce bra vector `x' to its canonical form
>     --
>     = compose cs bs
>       where
>           bs = basis x
>           cs = [x <> toKet b | b <- bs]

> braBasis :: Ord a => Bra a -> [Bra a]
>     --
>     -- List of unique basis of the bra vector
>     --
> braBasis BraZero        = []
> braBasis (Bra b)        = [Bra b]
> braBasis (_ :<| x)     = [x]
> braBasis (b1 :<+ b2)   = nub (braBasis b1 ++ braBasis b2)

> toKet :: Ord a => Bra a -> Ket a
>     --
>     -- Convert from bra to ket vector
>     --
> toKet (Bra k)            = Ket k
> toKet (x :<+ y)        = toKet x :+> toKet y
> toKet (p :<| Bra k)    = (conjugate p) :|> Ket k

> instance Ord a => DiracVector (Bra a)  where
>     zero          = BraZero
>     add           = (<+)
>     scale         = (<|)
>     reduce        = reduceBra
>     basis         = braBasis
>     components x  = [x <> toKet e | e <- basis x]
>     compose xs ks = foldl1 (:<+) [fst z :<| snd z  | z <- zip xs ks]

>     norm BraZero  = 0
>     norm x        = sqrt $ realPart (x <> toKet x)


Bra and Ket spaces as inner product spaces
==========================================

Definition. A complex vector space V is an inner product space if with
every pair of elements x ,y from V there is associated a unique inner
(or scalar) product \< x | y \> from C, such that

9:  < x | y >          = < y | x >*
10: < a x | b y >      = a* b < x | y >
11: < z | a x + b y >  = a < z | x > + b < z, y >
where
a, b, c are the complex scalars

We submit that the dual ket and bra spaces are inner product spaces,
providing that the inner product is defined by the operator \<\> given
below:

> (<>) :: Ord a => Bra a -> Ket a -> Scalar
>     --
>     -- Inner product, or the "bra-ket" product
>     --
> BraZero       <> _              = 0
> _             <> KetZero        = 0
> Bra i         <> Ket j          = d i j
> (p :<| x)     <> (q :|> y)      = p * q * (x <> y)
> (p :<| x)     <> y              = p * (x <> y)
> x             <> (q :|> y)      = q * (x <> y)
> x             <> (y1 :+> y2)    = (x  <> y1) + (x <> y2)
> (x1 :<+ x2)   <> y              = (x1 <> y)  + (x2 <> y)

> d :: Eq a => a -> a -> Scalar
> d i j
>     --
>     -- Classical Kronecker's delta
>     -- for instances of Eq class
>     --
>     | i == j    = 1
>     | otherwise = 0

The expressions below illustrate the definitions 9-11. They are all
true.

9:  (toBra x <> y) == conjugate (toBra y <> x)
10: (toBra (a |> x) <> (b |> y)) == (conjugate a)*b*(toBra x <> y)
11: (toBra z <> (a |> x +> b |> y)) == a*(toBra z <> x) + b*(toBra z <> y)
where
x = (2 :+ 3) |> Ket 2
y = ((1:+2) |> Ket 3) +> Ket 2
z = Ket 2 +> Ket 3
a = 2:+1
b = 1


Linear operators
================

Linear operators, or simply operators, are functions from vector in
representation a *a* to vector in representation *b*

a :: Ket a -> Ket b

although quite often the operations are performed on the same representation.
The linear operators A are defined by

A (c1 | x > + c2 | y > ) = c1 A | x > + c2 A | y >

We will describe variety of special types of operators, such as inverse,
unitary, adjoint and hermitian. This is not an accident that the names of those
operators resemble names from matrix calculus, since Dirac vectors and operators
can be viewed as matrices.

With the exception of variety of examples, no significant amount of Haskell code
will be added here. This section is devoted mainly to documentation; we feel
that it is important to provide clear definitions of the operators, as seen from
the Haskell perspective.  Being a strongly typed language, Haskell might not
allow for certain relations often shown in traditional matrix calculus, such as

A = B

since the two operators might have in fact two distinct signatures. In matrix
calculus one only compares tables of unnamed numbers, while in our Haskell
formalism we compare typed entieties. For this reason, we will be threading
quite slowly here, from one definition to another to assure that they are
correct from the perspective of typing rules of Haskell.

Operator notation
-----------------

The notation

| y > = A | x >

is pretty obvious: operator A acting on vector | x \> produces vector | y \>. It
is not obvious though whether both vectors use the same representation. The
Haskell version of the above clarifies this point, as in this example:

y = a >< x
where
a :: Ket Int -> Ket (Int, Int)
a = ......

In this case it is seen the two vectors have distinct representations.  The
operator \>\< will be explained soon but for now treat is as an application of
an operator to a vector, or some kind of a product of the two.

The above can be also written as

| y > = | A x >

where the right hand side is just a defining label saying that the resulting
vector has been produced by operator A acting on | x \>.

Linear operators can also act on the bra vectors

< y | = < x | A
<---

providing that they have correct signatures. This postfix notation though is a
bit awkward, and not supported by Haskell. To avoid confusion we will be using
the following notation instead:

< y | = < A x |

which says that bra y is obtained from ket y, where | y \> = | A x \>, as
before. In Haskell we will write it as

y = toBra $ a >< x

Renaming the representation
---------------------------

One simple example of an operator is *label "new"* which renames a vector
representation by adding extra label *"new"* in the basis vectors *Ket a*. Silly
as it sounds, this and other similar re-labeling operations can be actually
quite useful; for example, we might wish to distinguish between old and new
bases, or just to satisfy the Haskell typechecker.


label :: (Ord a, Ord b) => b -> Ket a -> Ket (b, a)
label i (Ket a) = Ket (i, a)
label i x       = (label i) >< x

Closure formula, or identity operator
-------------------------------------

Although the general Dirac formalism often refers to abstract vectors | x \>,
our implementation must be more concrete than that -- we always represent the
abstract vectors in some basis of our choice, as in:

| x > = ck | k >   (sum over k)

To recover the component c~k~ we form the inner product

ck = < k | x >

Putting it back to the previous equation:

| x > = < k | x > | k >      (sum over k)
= | k > < k | x >
= Id | x >
where
Id = | k > < k |        (sum over k)

we can see that the vector | x \> has been abstracted away. The formula says
that vector | x \> can be decomposed in any basis by applying identity operator
Id to it. This is also known as a closure formula. Well, Haskell has the "id"
function too, and we could apply it to any ket, as in:

id (Ket 1 +> 10 |> Ket 2) ==> | 1 > + 10 | 2 >

but Haskell's "id" does not know anything about representations; it just gives
us back the same vector | x \> in our original representation.

We need something more accurately depicting the closure formula | k \> \< k |,
that would allow us to change the representation if we wanted to, or leave it
alone otherwise. Here is the *closure* function and coresponding operator (\>\<)
that implement the closure formula for a given *operator*.


> closure :: (DiracVector a, DiracVector b) => (a -> b) -> a -> b
> closure operator x =
>    compose' (components x) (map operator (basis x))
>      where
>         compose' xs ks
>               | length xs == 0 = zero
>               | otherwise = foldl1 add (zipWith scale xs ks)

> (><) :: (DiracVector b, DiracVector a) => (a -> b) -> a -> b
> operator >< x = closure operator x

Changing the representation
---------------------------

The silly *label* function found in the comment of the section 6.1 uses in fact
the closure relation. But we could define is simpler than that:


> label :: t -> Ket t1 -> Ket (t, t1)
> label i (Ket x) = Ket (i, x)

and then apply a closure to a vector x, as in:

closure (label 0) (Ket 2 +> 7 |> Ket 3)
==> 1.0 |(0,2)> + 7.0 |(0,3)>

Somewhat more realistic example involves "rotation" of the old basis with
simulaneous base renaming:


> rot :: Ket Int -> Ket (Int, Int)
> rot (Ket 1) = normalize $ Ket (1,1) +> Ket (1,2)
> rot (Ket 2) = normalize $ Ket (1,1) +> (-1) |> Ket (1,2)
> rot (Ket _) = error "exceeded space dimension"

The example function *rot* assumes transformation from two-dimensional basis
[| 1 \>, | 2 \>] to another two-dimensional basis [| (1,1) \>, | (1,2) \>] by
expressing the old basis by the new one. Given this transformation we can apply
the closure to any vector | x \> represented in the old basis; as a result we
will get the same vector | x \> but represented in the new basis.

rot >< (Ket 1 +> 7 |> Ket 2) ==>
5.65685 |(1,1)> + -4.24264 |(1,2)>

Implementation of the operator equation A | x \> = | y \>
---------------------------------------------------------

The Haskell implementation of the closure formula is not just a useless
simulation of the theoretical closure - it is one of the workhorses of the
apparatus employed here.

We will be using linear operators to evaluate equations like this:

| y > = A | x >

The resulting vector | y \> can have either the same representation as | x \> or
different - depending on the nature of operator A. The most general type of A is

Ket a -> Ket b

but more often than not the basis will be the same as before. But how we define
the operator A itself? The best way is to specify how it acts on the base
vectors | k \>. If we can chose as our basis the eigenvectors of A this would be
even better, because the definition of A would be then extremely simple. After
inserting the identity | k \>\< k | between the operator A and vector | x \> in
the above equation one gets

| y > = A | k > < k | x >            (sum over k)

This will be implemented in Haskell as:

y = a >< x

The closure formula will take care of the rest and it will produce the result
| y \> . The examples previously given do just that. One caveat though: since
operator A will only be defined for the basis, but not for other vectors,
skipping the closure formula and coding directly

y = a' x

is not advisable. This will certainly fail for vectors other than
basis unless one makes extra provisions for that. This is what we did
in module Momenta, before we had the closure support ready. Using the
closure is safe and this is the way to go!

Inverse operator
----------------

An operator B = A^-1^ that inverses the equation

| y > = A | x >
y   = a >< x -- where a :: Ket a -> Ket b

into

| x > = B | y >
x   = b >< y -- where b :: Ket b -> Ket a

is called the inverse operator.

For example, the inverse operator to the operator *label i* is:


> label' :: (Ord a, Ord b) => Ket (a, b) -> Ket b
> label' (Ket (_, x)) = Ket x

It is easy to check that applying the operator A and its inverse A^-1^ in
succession to any ket | x \> one should obtain the same vector | x \> again, as
in:

A-1 A | x > = | x >

```haskell
label' >< (label 0 >< x) == x
    where
        x = Ket 1 +> 10 |> Ket 7

==> True
```

Once again, notice the omnipresent closure operator in Haskell implementation.
Tempting as it might be to implement the above example as

```bad-haskell
(label' . label 0) >< x == x
    where
        x = Ket 1 +> 10 |> Ket 7

==> True
```

this is not a recommended way. Although this example would work, but a similar
example for *rotation* operations would fail in a spectacular way. The correct
way is to insert the closure operator between two rotations:

```haskell
rot' >< (rot >< x) == x
    where
        x = Ket 1 +> 10 |> Ket 2

==> True
```

where the inverse operator *rot'* is defined below:


> rot' :: Ket (Int, Int) -> Ket (Int)
> rot' (Ket (1,1)) = normalize $ Ket 1 +> Ket 2
> rot' (Ket (1,2)) = normalize $ Ket 1 +> (-1) |> Ket 2
> rot' (Ket (_,_)) = error "exceeded space dimension"

Matrix representation of an operator
------------------------------------

The scalar products

< k | A l' > = < k | A | l' >

such that | k \> and | l' \> are the base vectors (in general belonging to two
different bases), form a transformation matrix Akl'.

In Haskell this matrix is formed as

```psuedo-haskell
k <> a >< l'
    where
        k  = ... :: Bra b
        l' = ... :: Ket a
        a  = ... :: Ket a -> Ket b
```

Adjoint operator
----------------

Our definition of adjoint operator is different than that in theory of
determinants. Many books, not necessarily quantum mechanical oriented, refer to
the latter as *classical adjoint operator*.

With every linear operator A we can associate an adjoint operator B = A^+^, also
known as Hermitian conjugate operator, such that equality of the two scalar
products

< A+ u | x > = < u | A x >

holds for every vector | u \> and | x \>. In Haskell notation the above can be
written as:

```psuedo-haskell
(toBra (b >< u) <> x) == toBra u <> a >< x
    where
        a = ... :: Ket a -> Ket b
        b = ... :: Ket b -> Ket a
        x = ... :: Ket a
        u = ... :: Ket b
```

For example, the operator *rot'* is adjoint to operator *rot*

```haskell
(toBra (rot' >< u) <> x) == (toBra u <> rot >< x)
    where
        x = Ket 1 +> 10 |> Ket 2
        u = Ket (1,1) +> 4 |> Ket (1,2)

==> True
```

It can be shown that

(A+)+ = A

Matrix A^+^ is conjugate transposed to A, as proven below

= A+kl'
= < k | A+ | l' >
= < k | A+ l' >
= < A+ l' | k >*
= < l' | A | k >*
= A*l'k

Unitary operator
----------------

Unitary transformations preserve norms of vectors. We say, that the norm of a
vector is invariant under unitary transformation. Operators describing such
transformations are called unitary operators.

< A x | A x > = < x | x >

The example of this is rotation transformation, which indeed preserves the norm
of any vector x, as shown in this Haskell example

```haskell
(toBra u <> u) == (toBra x <> x)
    where
        u = rot >< x
        x = Ket 1 +> 10 |> Ket 2

==> True
```

Inverse and adjoint operators of unitary operators are equal

A-1 = A+

which indeed is true for our example operator *rot*.

Computation of the adjont operators A^+^ from A is quite easy since the process
is rather mechanical, as described in the previous section. On the other hand,
finding inverse operators is not that easy, with the exception of some simple
cases, such as our example 2D rotation. It is therefore important to know
whether a given operator is unitary, as this would allow us to replace inverse
operators by adjoint operators.

Hermitian operator
------------------

A Hermitian operator is a self adjoint operator; that is

< A u | x > = < u | A x >

Another words: A^+^ = A.

Notice however, that this relation holds only for the vectors in the same
representation, since in general the operators A and A^+^ have distinct
signatures, unless types a, b are the same:

```haskell
a  :: Ket a -> Ket b -- operator A
a' :: Ket b -> Ket a -- operator A+
```

Elements of hermitian matrices must therefore satisfy:

Aij = (Aji)*

In particular, their diagonal elements must be real.

Our example operator *rot* is not hermitian, since it describes transformation
from one basis to another. But here is a simple example of a hermitian operator,
which multiplies any ket by scalar 4. It satisfies our definition:

```haskell
(toBra (a >< u) <> x) == (toBra u <> a >< x)
    where
        a v = 4 |> v

x = Ket 1 +> Ket 2
u = Ket 2

==> True
```

Here is a short quote from [3].

Why do we care whether an operator is Hermitian? It's because of a few theorems:

-   1. The eigenvalues of Hermitian operators are always real.
-   2. The expectation values of Hermitian operators are always real.
-   3. The eigenvectors of Hermitian operators span the Hilbert space.
-   4. The eigenvectors of Hermitian operators belonging to distinct eigenvalues
    are orthogonal.

In quantum mechanics, these characteristics are essential if you want to
represent measurements with operators. Operators must be Hermitian so that
observables are real. And, you must be able to expand in the eigenfunctions -
the expansion coefficients give you probabilities!


Showing kets and bras
=====================

Lastly, here are show functions for pretty printing of Dirac vectors.


> instance (Show a, Eq a, Ord a) => Show (Ket a)  where
>     showsPrec _ KetZero   = showString "| Zero >"
>     showsPrec n (Ket j)   = showString "|" . showsPrec n j . showString ">"
>     showsPrec n (x :|> k) = showsScalar n x . showsPrec n k
>     showsPrec n (j :+> k) = showsPrec n j . showString " + " . showsPrec n k

> instance (Show a, Eq a, Ord a) => Show (Bra a)  where
>     showsPrec _ BraZero   = showString "< Zero |"
>     showsPrec n (Bra j)   = showString "<" . showsPrec n j . showString "|"
>     showsPrec n (x :<| k) = showsScalar n x . showsPrec n k
>     showsPrec n (j :<+ k) = showsPrec n j . showString " + " . showsPrec n k


> showsScalar :: (Show t, RealFloat t) => Int -> Complex t -> String -> String
> showsScalar n x@(a :+ b)
>     | b == 0    = showString " " . showsPrec n a
>     | otherwise = showString " (" .showsPrec n x . showString ")"


Data Tuple for tensor products
==============================

A state vector of several subsystems is modelled as a ket parametrized by a type
variable Tuple, which is similar to ordinary () but is shown differently. Tensor
product of several simple states leads to deeply entangled structure, with many
parenthesis obstructing readability.  What we really want is a simple notation
for easy visualization of products of several states, as in:

Ket 1 *> Ket (2, 1) * Ket '+' ==> | 1; (2,1); '+' >

See module Momenta for practical example of tensor products of vector spaces.


> data Tuple a b =  a :* b
>     deriving (Eq, Ord)

> instance (Show a, Show b) => Show (Tuple a b) where
>     showsPrec n (a :* b) = showsPrec n a . showString "; " . showsPrec n b


Messing Around
==============

> data Operator k b =
>   OpIdent
>   | Ket k :>-< Bra b
>   | Operator k b :>+< Operator k b

> (>-<) :: Ket a -> Bra b -> Operator a b
> kt >-< br = kt :>-< br

> (>+<) :: (Ord k, Ord b, Eq k, Eq b)
>   => Operator k b -> Operator k b -> Operator k b
> op1 >+< op2 = op1 :>+< op2
> -- If correct, we could do this (possible speed up evaluation)
> -- (k1 :>-< b1) >+< (k2 :>-< b2) = (k1 +> k2) >-< (b1 <+ b2)

> (>*<) :: (Ord k1, Ord k2, Ord b1, Ord b2, Eq k1, Eq k2, Eq b1, Eq b2)
>   => Operator k1 b1 -> Operator k2 b2 -> Operator (Tuple k1 k2) (Tuple b1 b2)
> (k1 :>-< b1) >*< (k2 :>-< b2) = (k1 *> k2) >-< (b1 <* b2)

> (>*) :: (Ord bi, Ord bf, Eq bi, Eq bf)
>   => Operator bf bi -> Ket bi -> Ket bf
> _              >* KetZero = KetZero
> --OpIdent        >* x       = 
> (op1 :>+< op2) >* x       = (op1 >* x) +> (op2 >* x)
> (kts :>-< brs) >* x       = (brs <> x) |> kts

> (*<) :: (Ord bi, Ord bf, Eq bi, Eq bf)
>   => Bra bi -> Operator bi bf -> Bra bf
> BraZero *< _              = BraZero
> --x       *< OpIdent        = x
> x       *< (op1 :>+< op2) = (x *< op1) <+ (x *< op2)
> x       *< (kts :>-< brs) = (x <> kts) <| brs

> instance (Show a, Show b, Ord a, Ord b, Eq a, Eq b)
>   => Show (Operator a b) where
>   showsPrec n OpIdent = showString "|I><I|"
>   showsPrec n (op1 :>+< op2) =
>       showsPrec n op1 . showString " >+< " . showsPrec n op2
>   showsPrec n (k :>-< b) = showsPrec n k . showsPrec n b


References
==========

-   [1] Jerzy Karczmarczuk, Scientific computation and functional
programming, Dept. of Computer Science, University of Caen,
France, Jan 20, 1999,
[http://www.info.unicaen.fr/\~karczma/](http://www.info.unicaen.fr/~karczma/)
-   [2] Jan Skibinski, Collection of Haskell modules, Numeric Quest
Inc.,
[http://www.numeric-quest.com/haskell/"](http://www.numeric-quest.com/haskell/)
-   [3] Steven Pollock, University of Colorado, [Quantum Mechanics,
Physics 3220 Fall 97, lecture
notes](http://www.colorado.edu/physics/phys3220/3220_fa97/notes/notes_table.html)


Copyright and license
=====================

--
-- Copyright:
--
--      (C) 2000 Numeric Quest, All rights reserved
--
--      Email: jans@numeric-quest.com
--
--      http://www.numeric-quest.com
--
-- License:
--
--      GNU General Public License, GPL
--
