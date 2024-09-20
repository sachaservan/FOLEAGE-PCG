# $\mathbb{F}_4$ OLEAGE PCG

A prototype implementation of the $\mathbb{F}_4$ OLEAGE Pseudorandom Correlation Generator (PCG) in C.
See the [paper](https://eprint.iacr.org/2024/429.pdf) for details.

## Organization

The [libs/](libs/) folder contains the implementation of the (parallel) FFT and [ternary DPF](https://github.com/sachaservan/tri-dpf) (a submodule) which are both used extensively in the PCG construction.

The [src/](src/) folder contains the code for PCG implementation.

- [src/test.c] closely follows Figure 1 from the [paper](https://eprint.iacr.org/2024/429.pdf) and uses parameters `c=4, t=27` which simplifies the implementation by having the number of blocks `t` be a power of 3.
- [src/bench.c] implements two benchmarks and can be used to reproduce Table 6 in the paper.

## Dependencies

These dependencies are required by the [ternary DPF](https://github.com/sachaservan/tri-dpf) submodule.

- OpenSSL
- GNU Make
- Cmake
- Clang

## Getting everything to run (tested on Ubuntu, CentOS, and MacOS)

| Install dependencies (Ubuntu):         | Install dependencies (CentOS):              |
| -------------------------------------- | ------------------------------------------- |
| `sudo apt-get install build-essential` | `sudo yum groupinstall 'Development Tools'` |
| `sudo apt-get install cmake`           | `sudo yum install cmake`                    |
| `sudo apt install libssl-dev`          | `sudo yum install openssl-devel`            |
| `sudo apt install clang`               | `sudo yum install clang`                    |

On MacOS, use [homebrew](https://brew.sh/) to install dependencies.
`cmake` and `clang` can be installed via `xcode-select --install`.
OpenSSL can be installed via `brew install openssl` or manually.

## Running tests and benchmarks

Test:

```
git submodule update --init --recursive
make
./bin/pcg --test
```

Benchmarks:

```
git submodule update --init --recursive
make
./bin/pcg --bench
```

DPF benchmarks:
See the [DPF repository](https://github.com/sachaservan/tri-dpf).

SPFSS benchmarks:
Since the SPFSS benchmarks are specific to FOLEAGE, we provide a special test file `spfss_test.c` which can be used to benchmark the DPF implementation. To do so, run:

```
cd libs
mv tri-dpf/src/test.c tri-dpf/src/test.old
cp spfss_test.c tri-dpf/src/test.c
cd tri-dpf
make && ./bin/test
```

FFT benchmarks

```
cd libs/fft
make && ./bin/fft
```

## Parameter Selection

The parameters `c` and `t` can be computed using the [SageMath parameter selection script](https://github.com/mbombar/estimator_folding) (also available as a submodule in `scripts/parameters_selection`).
We provide reasonable choices of `c` and `t` in Table 2 of [the paper](https://eprint.iacr.org/2024/429.pdf).
In particular, our benchmarks use `(c=4, t=27)` as a conservative parameter choice and `(c=3,t=27)` as an aggressive parameter choice, when targeting at least $\lambda=128$ bits of security.

## Future development

The current prototype implementation can be extended in several ways.
TODOs are left in-line, however, the broad strokes include:

- [ ] Unit tests for the FFT (currently only checked by hand on a small instance).
- [ ] Modularize the PCG construction and tests (currently the test is one monolithic block).
- [ ] More efficient SIMD-based implementation of the FFT packing (a matrix transpose) which currently is implemented using the naive approach and results in the computational performance bottleneck.

## Citation

```
@inproceedings{foleage,
  author       = {Maxime Bombar and
                  Dung Bui and
                  Geoffroy Couteau and
                  Alain Couvreur and
                  Clément Ducros and
                  Sacha Servan-Schreiber},
  title        = {FOLEAGE: $\mathbb{F}_{4}$OLE-Based Multi-Party
                  Computation for Boolean Circuits},
  note         = {\url{https://eprint.iacr.org/2024/429}},
  editor       = {Kai-Min Chung and
                  Yu Sasaki},
  booktitle    = {Advances in Cryptology - {ASIACRYPT} 2024 - 30th
                  International Conference on the Theory and
                  Application of Cryptology and Information Security,
                  Kolkata, India, December 9-13, 2024 %
                  },
  publisher    = {Springer},
  year         = {2024},
}
```

## ⚠️ Important Warning

<b>This implementation is intended for _research purposes only_. The code has NOT been vetted by security experts.
As such, no portion of the code should be used in any real-world or production setting!</b>
