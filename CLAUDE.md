# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
# Linux (primary)
make

# Haiku OS
make -f Makefile.haiku
```

The Linux build uses `gcc -O3 -march=native -std=gnu89 -ffast-math -fopenmp` and links against `-lm -lfftw3`. The binary is `./gabord`.

**Performance:** On Haiku, set `OMP_NUM_THREADS=4` before running, or call `omp_set_num_threads(4)` early in `main()`.

## Run

```bash
./gabord
```

Reads `data.bin` (raw `double` array) from the current directory and writes decomposition results to stdout. (FFTW now uses `FFTW_ESTIMATE`, so there is no longer a `fftw_wisdom.dat` file — see *FFTW integration* below.)

## Architecture

This is a **Gabor Matching Pursuit** signal decomposition library, originally based on Mallat/Zhang (1993), adapted for continuous EEG processing at Johns Hopkins.

### Core algorithm flow

```
main()
  → gabord()            # top-level entry point (gabord.c)
      → GaborBuildFilter()   # constructs Gabor filter bank (gb_filter.c)
      → GaborDecomp()        # Gabor transform via FFTW (gb_decomp.c)
      → GaborBuildBook()     # one iteration of Matching Pursuit (gb_buildbook.c)
          → GaborGetMaxFrmTrans()   # find best atom across all octaves (gb_oper.c)
          → GaborGetResidue()       # subtract selected atom
          → UpdateGabor()           # update transform coefficients
          → GaborUpdateFourier()    # update stored Fourier transform
```

`gabord()` is designed to be called repeatedly (e.g., from Matlab via mex, or per sliding window). It uses a `gabord_initialised` guard (file-scope, threadprivate) to allocate filter banks and signal arrays only on the first call per thread; subsequent calls reuse filters but reset the book. See *OpenMP parallelism & thread safety* for concurrent use.

### Data structures (mpp.h / gabsignals.h)

- **`GABSIGNAL`** — pointer to `struct gabsignal`. Stores a split-complex array: `values[0..size/2-1]` = real part, `values[size/2..size-1]` = imaginary part. Many routines depend on this layout.
- **`FILTER`** — typedef alias for `GABSIGNAL`; filter bank entries share the same struct.
- **`BOOK`** — linked list of `WORD`s, each holding an atom (octave, frequency id, position, coefficient, phase).
- **`cur_book`**, **`cur_filter`**, **`cur_transform`** — macros into `library[]`, `filter[]`, `transform[]` arrays indexed by `Current_Book` (always 0 in current usage).

### Key global state (gabord.c)

Most algorithm parameters live as file-scope globals: `cur_MinOctave`, `cur_MaxOctave`, `cur_SOT/SOF`, `pfG/pfC/pfB/pfCE1/pfCE2/pfCE3`, `pnAep/pnIep`. These are lazily allocated and freed when signal size changes.

### FFTW integration (gb_decomp.c)

Plans are cached in a `plan_cache[]` array (up to 32 entries) and created with
`FFTW_ESTIMATE`. The transforms are tiny (~1.4% of runtime), so measured plans bought
nothing while `FFTW_MEASURE` had two drawbacks for this workload: it is **non-deterministic**
(picks codelets by wall-clock benchmarking, so results vary run-to-run at the 1e-13 level)
and it relied on a CWD-shared `fftw_wisdom.dat` file that races across concurrent processes.
`FFTW_ESTIMATE` is deterministic and needs no wisdom file. `GaborFFTCleanup()` is registered
via `atexit()`.

`plan_cache`/`plan_cache_count` are **threadprivate** (each thread owns its scratch buffers,
so concurrent `fftw_execute` calls never collide), and the FFTW planner — which is **not**
thread-safe — is serialised with `#pragma omp critical(fftw_planner)`. Execution stays fully
parallel.

### OpenMP parallelism & thread safety

`gabord()` is safe to call **concurrently from multiple threads** (e.g., one channel per
thread). All mutable per-call state is `threadprivate` — see the `#pragma omp threadprivate`
lists in `mpp.h` and `gabord.c`, plus `plan_cache` (gb_decomp.c), `indexG1/indexG2`
(gb_oper.c) and the cache vars `nL/maxlh` (gabord.c). **Each worker thread must call
`gabord_reset()` once before its first `gabord()` call** (threadprivate state is otherwise
undefined in new threads); `gabord_cleanup()` frees a thread's state.

Two parallelism modes, and they are mutually exclusive by design:

- **Internal** — five `#pragma omp parallel for` sites (`GaborGetMaxFrmTrans` in gb_oper.c;
  `GaborUpdateFourier` ×2 in gb_buildbook.c; `UpdateGabor` ×2 in update.c) parallelise a
  single decomposition over octaves / frequency bins. This helps only a lone call: it covers
  ~30% of the work (the `genericLoop` matching-pursuit update in update.c is serial), so its
  Amdahl ceiling is ~1.4×, and on small windows the overhead is barely worth it.
- **External (preferred for batch)** — call `gabord()` on many independent signals from an
  outer `omp parallel` region. Each of the five internal sites is guarded with
  `if(!omp_in_parallel())`, so they **auto-disable** under outer parallelism, avoiding nested
  teams. This is how `maf_gabord --jobs N` scales across channels.

### Signal size constraint

Input size is rounded up to the next power of 2 inside `gabord()`. The maximum supported size is `SIG_SIZE = 32768` (defined in `gabsignals.h`).

### Stop criteria (gabord.c)

Controlled by four globals (set before calling `gabord()`):
- `citer` / `max_num_iter` — stop after N atoms
- `cpct` / `th_pct` — stop when book energy reaches a % of signal energy
- `cth` / `thatomnrj` — stop when last atom energy falls below threshold
- `ccoh` — coherence-based criterion (only valid when `cur_SOT==cur_SOF==1`)

If none are set, `citer=1, max_num_iter=1000` is the default.
