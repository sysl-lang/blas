# blas

The Basic Linear Algebra Subprograms for sysl, through the CBLAS interface: all 146 routines of the
reference CBLAS at `s`, `d`, `c` and `z`, and over them column-major matrix views and strided vector
views that check every shape before BLAS sees the data.

*An independent sysl package for BLAS. It is not affiliated with or endorsed by the authors of the
reference BLAS (Netlib), of OpenBLAS, or of Apple's Accelerate, or by their trademark holders.*

```sysl
import sh.sysl.blas.{MatrixView, VectorView, BlasError, NoTrans, Trans, gemm, gemv, dot}

main() -> Result[unit, BlasError]
    var a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]        // 2 × 3, column-major: [[1, 3, 5], [2, 4, 6]]
    var c = [0.0, 0.0, 0.0, 0.0]

    val m = MatrixView.of(a[..], 2, 3)?

    gemm(NoTrans, Trans, 1.0, m, m, 0.0, MatrixView.of(c[..], 2, 2)?)?
    print(c)                                       // A·Aᵀ: [35, 44, 44, 56]

    var x = [1.0, 1.0, 1.0]
    var y = [0.0, 0.0]

    gemv(NoTrans, 1.0, m, VectorView.of(x[..])?, 0.0, VectorView.of(y[..])?)?
    print(y)                                       // the row sums: [9, 12]

    print(dot(m.row(0), VectorView.strided(x[..], 3, -1)?)?)   // 1 + 3 + 5 = 9
    Ok(())
```

```
dependencies {
  blas { git = "github.com/sysl-lang/blas", version = "0.1.0" }
}
```

The coordinate names the package, `blas`; the module a program imports is `sh.sysl.blas`.

## Which BLAS

**Apple's Accelerate on macOS, OpenBLAS everywhere else** — the choice nex's compiler makes when it
links a program whose `@` goes to BLAS, and the same LP64 `cblas_*` symbols it calls. A product
computed through this package and the same product in a compiled nex program run the same routine
in the same library, so their results are bit-identical.

| | macOS | elsewhere |
|---|---|---|
| default (`platform`) | Accelerate: `@link("cblas")`, which the SDK's `libcblas.tbd` resolves to `Accelerate.framework`'s `libBLAS.dylib`; nothing to install | `-lopenblas` from the system's search path (`apt install libopenblas-dev`) |
| `features = [openblas]` | Homebrew's OpenBLAS, found through `openblas.pc` | OpenBLAS, found through `openblas.pc` |

Homebrew's OpenBLAS is keg-only, so `pkg-config` finds it only with
`PKG_CONFIG_PATH=/opt/homebrew/opt/openblas/lib/pkgconfig`. The two features are a `choices` group
(`library`), so exactly one is on.

```
sysl test .
PKG_CONFIG_PATH=/opt/homebrew/opt/openblas/lib/pkgconfig sysl test . --features openblas
```

## What it is

Two layers. **`sh.sysl.blas.c`** is CBLAS one declaration per routine — pointers, `i32` extents and
increments, the CBLAS enumerations as `i32` constants (`COL_MAJOR`, `TRANS`, `UPPER`, …) — and checks
nothing. **`sh.sysl.blas`** is what a program imports:

| | |
|---|---|
| **`MatrixView[T]`** | a column-major matrix over storage it does not own: `ptr`, `rows`, `cols`, `ld` (element `(i, j)` at `ptr + i + j·ld`). `of(xs, rows, cols)` (dense), `strided(xs, rows, cols, ld)`, `from_ptr(p, rows, cols, ld)`; `block(r, c, rows, cols)` keeps the parent's `ld`; `row(i)`, `col(j)`, `at`, `set`, `place` |
| **`VectorView[T]`** | a strided vector: `ptr`, `len`, `inc`. `of(xs)`, `strided(xs, len, inc)`, `from_ptr(p, len, inc)`; `at`, `set`, `place`. **A negative increment follows the reference BLAS rule**: `ptr` is the lowest-addressed element and the vector runs backwards from the far end |
| **`Op`** | `NoTrans`, `Trans`, `ConjTrans` — how an operand enters a product |
| `dot`, `dotc` | `Σ xᵢyᵢ` (unconjugated: nex's `x @ y`) and `Σ conj(xᵢ)yᵢ` |
| `axpy`, `scal` | `y ← αx + y`, `x ← αx` |
| `nrm2`, `asum` | at the element's real width (`T::Real`); `asum` of a complex vector is BLAS's `Σ |re| + |im|` |
| `gemv`, `gemm` | `y ← α·op(A)·x + β·y`, `C ← α·op(A)·op(B) + β·C` |
| **`BlasError`** | `ShapeMismatch`, `ZeroIncrement`, `LeadingDimension`, `OutOfStorage`, `TooLarge` |
| `Element` | the dispatch bound: `f32`, `real`, `Complex[f32]`, `Complex[real]` |

**Every check happens before C sees the data.** A view is checked when it is made — its farthest
element inside the slice, `ld ≥ max(1, rows)`, a non-zero increment, every extent within CBLAS's
32-bit `int` — and a routine checks only that its operands' extents agree. So no call can reach the
library's own parameter error (Accelerate prints it and aborts).

**From an array that is not a sysl slice** — nex's arrays, a C library's — `from_ptr` makes a view of
the address of the first element with its extents and stride, without a copy; only the leading
dimension, the increment and the 32-bit range are checked there.

**How nex's `@` maps onto it.** `A @ B` is `gemm(NoTrans, NoTrans, 1, A, B, 0, C)` with `C` dense
(`ld = max(1, m)`); `transpose(A)` and `conj(transpose(A))` are `A` with `Trans` and `ConjTrans`; a
matrix and a vector is `gemv`, a vector and a matrix is `gemv` with the flag flipped, and
`x @ conj(transpose(Y))` is `gemm` with `x` as a `1 × k` matrix (`MatrixView.of(x, 1, k)`). Two
vectors are `dot`.

**An output must not share elements with an input**, which is BLAS's own rule. It is not checked:
two blocks of one matrix interleave in memory without sharing an element, and an address-range test
would refuse them.

## Tested

`sysl test .` — 59 tests, against Accelerate and, with `--features openblas`, against OpenBLAS 0.3.34.
The oracle is never BLAS read back against itself: hand-computed small cases on integer-valued
operands (exact, so compared with `==`), identities computed by separate calls (`A·I = A`,
`(AB)ᵀ = BᵀAᵀ` over every pair of transpose flags, `trsm` undoing `trmm` over all sixteen
side/triangle/transpose/diagonal combinations, `trsv` undoing `trmv`), and numpy's answers for a
9 × 7 by 7 × 8 product and a complex `Aᴴ·B`, on operands chosen so every partial sum is exact. Strided
vectors at increment 2 and at negative increments, sub-block views with a leading dimension larger
than the row count, every `Op`, and every refusal. Every one of the 146 routines is called once, which
proves each symbol resolves in the library linked.

Run under AddressSanitizer (`SYSL_EXTRA_CFLAGS="-fsanitize=address -g" sysl test .`): the sysl half
is instrumented (`instrumented: address (1 binary)`), the library is not.

## What is deliberately absent

- **Each library's own extensions** — Accelerate's `_spi`/`_singlecore` forms, `catlas_*`,
  `appleblas_*`; OpenBLAS's `axpby`, `omatcopy`, `gemmt`, `i?amin`, the bfloat16 routines and the
  rest. A program linked against the other library would not find them. `SURVEY.md` lists both.
- **The ILP64 interface** (64-bit `int`). Neither library exports it by default, and nex calls LP64.
- **Row-major views.** The raw layer takes `ROW_MAJOR`; the views are column-major because that is
  BLAS's, LAPACK's and nex's layout. A row-major matrix is its transpose, column-major.
- **Pleasant forms of the rest of level 2 and 3** (`symv`, `trsv`, `syrk`, `trsm`, …): reachable
  through `sh.sysl.blas.c` today.
