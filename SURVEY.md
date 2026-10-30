# Survey of CBLAS (Accelerate and OpenBLAS 0.3.34)

Every routine of the CBLAS interface, classified before the binding was written: **expressible**
(declares directly), **shim-able** (needs C of its own), or **needs the language to move**.

## Which routines

The two libraries export different sets, so the binding declares their intersection — the reference
(Netlib) CBLAS — and nothing a program linked against the other library would fail to find:

| | count | how it was listed |
|---|---|---|
| Accelerate (`libcblas.tbd` in the macOS SDK → `vecLib.framework/libBLAS.dylib`) | 168 `cblas_*` | `grep -o "_cblas_[a-z0-9_]*"` of the `.tbd` |
| OpenBLAS 0.3.34 (`libopenblas.dylib`) | 213 `cblas_*` | `nm -gU` |
| both | 147 | `comm -12` |
| **declared** | **146** | both, less `cblas_xerbla` |

Accelerate-only (21, not declared): the `_spi` / `_spi_mem_requirement` forms of `?gemm` and `?trsm`,
the `_singlecore` forms of `?gemm` and `?gemv`, `cblas_errprn`; also `catlas_*` and `appleblas_*`,
which `cblas.h` declares under other prefixes. OpenBLAS-only (66, not declared): `?axpby`, `?axpyc`,
`?omatcopy`/`?imatcopy`, `?geadd`, `?gemmt`/`?gemmtr`, `?gemm_batch`, `i?amin`/`i?max`/`i?min`,
`?amax`/`?amin`, `?sum`, the bfloat16 routines (`sbgemm`, `sbgemv`, `sbdot`, the converters), and
`cdotu`/`cdotc`/`zdotu`/`zdotc` returning a complex by value. `cblas_xerbla` and `cblas_errprn` are
the libraries' error reporters, called by them rather than by a program.

## Result

| bucket | count | what |
|---|---|---|
| expressible | 146 routines, 5 enumerations | everything |
| shim-able | 0 | — |
| needs the language to move | 0 | — |

No shim is needed and no sysl gap was found.

## The routines, by level

| level | routines (each at `s`, `d`, `c`, `z` unless a prefix is named) | count |
|---|---|---|
| 1 | `dot` (`s`, `d`), `sdsdot`, `dsdot`, `dotu_sub`, `dotc_sub` (`c`, `z`), `nrm2` (`s`, `d`, `sc`, `dz`), `asum` (`s`, `d`, `sc`, `dz`), `i?amax`, `swap`, `copy`, `axpy`, `scal`, `csscal`, `zdscal`, `rotg`, `rot` (`s`, `d`), `csrot`, `zdrot`, `rotm`, `rotmg` (`s`, `d`) | 50 |
| 2 | `gemv`, `gbmv`, `trmv`, `tbmv`, `tpmv`, `trsv`, `tbsv`, `tpsv` (all four); `symv`, `sbmv`, `spmv`, `ger`, `syr`, `spr`, `syr2`, `spr2` (`s`, `d`); `hemv`, `hbmv`, `hpmv`, `geru`, `gerc`, `her`, `hpr`, `her2`, `hpr2` (`c`, `z`) | 66 |
| 3 | `gemm`, `symm`, `syrk`, `syr2k`, `trmm`, `trsm` (all four); `hemm`, `herk`, `her2k` (`c`, `z`) | 30 |
| | | **146** |

## The shapes, and how each declares

| shape | where | how |
|---|---|---|
| `int` dimension, increment, leading dimension | every routine | `i32` — LP64, the interface both libraries export by default and the one nex's compiler calls |
| `enum CBLAS_ORDER` / `_TRANSPOSE` / `_UPLO` / `_DIAG` / `_SIDE` | levels 2 and 3 | `i32`, the values written as `val`s (101–102, 111–113, 121–122, 131–132, 141–142: the CBLAS standard's, identical in both headers) |
| real array `float *`, `double *` | everywhere | `*f32`, `*f64` |
| complex array or scalar `void *` | every `c`/`z` routine | `*Complex[f32]`, `*Complex[f64]` — two fields of the width, C's `float _Complex` layout |
| complex result | `?dotu_sub`, `?dotc_sub` | the out-parameter, a `*Complex[F]` to a local |
| real scalar by value | `alpha`, `beta`, `c`, `s` | `f32` / `f64` |
| `CBLAS_INDEX` result | `i?amax` | `i32` — Accelerate's is `int`, OpenBLAS's `size_t`; the index is in the low 32 bits of the return register in both |

## What forced a decision, though not a gap

- **No header is compiled.** The headers differ (Accelerate's `cblas.h` is under the SDK's
  `vecLib.framework`, deprecated for `ACCELERATE_NEW_LAPACK`; OpenBLAS's is keg-only under Homebrew,
  and `openblas.pc`'s `Cflags` carry `-Xpreprocessor -fopenmp`). The symbols are the same, so
  `c/c.sysl` declares them and the package compiles no C. A wrong declaration would show as a wrong
  answer, and the tests check every routine family against hand-computed values; every one of the 146
  is called once (`tests_symbols.sysl`), which proves each resolves in the library linked.
- **Accelerate is linked as `-lcblas`, not `-framework Accelerate`.** A sysl `@link` names a library,
  never a flag; the SDK's `libcblas.tbd` points at the same `libBLAS.dylib` inside
  `Accelerate.framework` that `-framework Accelerate` reaches (`otool -L` shows it). So the code that
  answers is the code nex's compiled programs call.
- **The legacy (non-`ACCELERATE_NEW_LAPACK`) symbols.** nex's compiler emits calls to the plain
  `cblas_?gemm` names, LP64; this binding calls the same symbols, so an interpreted product and a
  compiled one run the same routine. The new-header ILP64 forms are not declared.
- **A complex element dispatches through its width.** sysl's `impl` for a generic type covers every
  instantiation (`traits.md`, "An `impl` covers a generic type as a whole"), so `Complex[f32]` and
  `Complex[real]` cannot each carry their own `impl Element`. `impl[F: ComplexWidth] Element for
  Complex[F]` is one implementation, and `ComplexWidth` — implemented by `f32` and `real` — says which
  of the `c` and `z` routines a width is. That is the language's stated design, not a gap.
