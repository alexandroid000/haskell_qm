---
title: Quantum Mechanics in Haskell
...


Building on the work of Jan Skibinski and Jerzy Karczmarczuk, we hope to
demonstrate how Quantum Systems may be intuitively built and simulated using the
Haskell functional programming language. Many features of Haskell make it useful
for scientific computation.

The file [QuantumVector.lhs](QuantumVector.lhs) sets up basic Dirac notation for
manipulating Quantum systems, which we will build on to demonstrate some simple
systems. The file [scientific_computation.pdf](scientific_computation.pdf) talks
about some other useful features of Haskell for programming scientific systems.

The goal of using Haskell to build scientific systems is readability and
usability. Haskell allows us to define the mathematics of the system we're
setting up (a Hilbert space for quantum mechanics), using the notation we like
(Dirac notation for quantum). This makes it easier to take our analytical
reasoning about a system and apply it directly to the programs we create.

Start by reading the [quantum vector](QuantumVector.lhs) file. Having the
Haskell compiler installed will help a lot, as you can load the file into the
interpreter (GHCi) and play around with the math a bit. With this file loaded,
we can perform basic quantum math:

-   Create 'Bra's and 'Ket's which represent our quantum states
-   Add, subtract, scale, and tensor 'Bra's and 'Ket's together
-   Create superposition states and normalize them
-   Define a basis
-   Dot-product (inner product) a 'Bra' and a 'Ket' together
-   Define operators and apply them to states
-   Project a state onto another basis

The math is a bit hard to read, because we can't use the normal symbols `+`,
`-`, and `*`. These definitions are in the file towards the top too.

-   `|>` scale a ket
-   `<|` scale a bra
-   `*>` tensor product of two kets
-   `<*` tensor product of two bras
-   `+>` sum of two kets
-   `<+` sum of two bras
-   `<>` inner product
-   `><` closure (apply operator)
-   `a :+ b` complex number `a + bi`

Below are some examples of how to use the Quantum Vector module. To use it,
install `ghc` and `ghci` on your system (the Glasgow Haskell Compiler and the
interpreter for the ghc). At your command prompt, type `ghci` to launch the
interpreter, then load the `QuantumVector.lhs` script by issuing the command
`:l QuantumVector.lhs` to the interpreter.

```haskell
> :l QuantumVector.lhs

--  Make some Kets, add them
QuantumVector> Ket 3
|3>
QuantumVector> 2.0 |> Ket 3
2.0 |3>
QuantumVector> 2.0 |> Ket 3 +> 1.0 |> Ket 1
2.0 |3> + 1.0 |1>

--  Dot product behaves as expected
QuantumVector> Bra 3 <> Ket 3
1.0 :+ 0.0
QuantumVector> Bra 3 <> Ket 2
0.0 :+ 0.0

--  Make a superposition of states
QuantumVector> compose [1,2,3,4] [Ket 1, Ket 2, Ket 3, Ket 4]
1.0 |1> + 2.0 |2> + 3.0 |3> + 4.0 |4>

--  Can also normalize a superposition
QuantumVector> normalize $ compose [1,2,3,4] [Ket 1, Ket 2, Ket 3, Ket 4]
0.18257418583505536 |1> + 0.3651483716701107 |2> + 0.5477225575051661 |3> + 0.7302967433402214 |4>

--  Labels can eb anything (as long as it's orderable), here we use Tuples
QuantumVector> Ket (2,3)
|(2,3)>
QuantumVector> Ket (2,3) +> Ket (2,2)
1.0 |(2,3)> + 1.0 |(2,2)>
QuantumVector> normalize $ Ket (2,3) +> Ket (2,2)
0.7071067811865475 |(2,3)> + 0.7071067811865475 |(2,2)>

--  Can extract a basis from Kets added together
QuantumVector> basis $ Ket 3
[|3>]
QuantumVector> Ket 3 +> Ket 2 +> Ket 5
1.0 |3> + 1.0 |2> + 1.0 |5>
QuantumVector> basis $ Ket 3 +> Ket 2 +> Ket 5
[|3>,|2>,|5>]
QuantumVector> basis $ Ket 3 +> Ket 2 +> 2.0 |> Ket 5
[|3>,|2>,|5>]

--  Can tensor together Kets to make a new bigger space
QuantumVector> Ket 3 *> Ket 2 *> Ket 5
|3; 2; 5>
```

One important part of the formalism in `QuantumVector` is the closure operator,
which allows us to apply operators to bras and kets, like so: `A | x \> = | y \>`

To implement this operation, we will define a general closure operator, `><`,
which works by inserting the identity (`| k \> \< k |`) between `A` and `x`,
where `k` is the basis of `x`. Look at the definition of the closure operator:

```haskell
closure :: (DiracVector a, DiracVector b) => (a -> b) -> a -> b
closure operator x =
   compose' (components x) (map operator (basis x))
     where
        compose' xs ks
              | length xs == 0 = zero
              | otherwise = foldl1 add (zipWith scale xs ks)

(><) :: (DiracVector b, DiracVector a) => (a -> b) -> a -> b
operator >< x = closure operator x
```

The first line is the type signature. Type signatures are optional in Haskell,
but they make it easier to understand what functions are supposed to do. Let's
break it down:

- `(DiracVector a, DiracVector b)` tells us that `a` and `b` in the type
  signature are bras or kets, which are described by the general typeclass
  `DiracVector`.
- After the `=>`, we describe the actual arguments and results of the function.
- `(a -> b)` tells us that the first argument to the closure function is a
  function, which transforms `a` into `b`. This is our operator.
- The second argument, after the next arrow, is just an `a`, so a bra or a ket.
  This is the vector our operator is applied to.
- The last part of the type signature, after the last `->`, is what the function
  `closure` returns. In our case, it's `b`, which is the type that the operator
  produces when it operates on and `a` type.

So this type signature describes a function which takes an operator `A`, a vector
`x`, and returns the result of `A` applied to `x`.

The nitty-gritty of how we implement this in Haskell is less important, but we
can see in the first line of the actual function definition that we separate out
the `basis` and `components`, apply the operator to the basis of `x`, and then
recompose the vector. The last two lines of the code block just redefine the
`closure` function to look pretty, in the Dirac notation. The [quantum
vector](QuantumVector.lhs) file defines a rotation operator and shows how to
apply it with the closure operator, so take a look at that section ("Changing
the representation").

For the talk, we will review this module and answer any questions about it, then
we will move on to representing the annihilation and creating operators in
Haskell. We'll start by making them for the quantum harmonic oscillator, then
move on towards generalizing them.
