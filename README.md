# Standalone FoXS with GAMB

Build the default CPU executable:

```sh
make
```

Run the regression test:

```sh
make test
```

## Optional CUDA acceleration

The CUDA build accelerates the weighted pair-distance distributions used for
normal and partial SAXS profiles. The pair kernel evaluates only the upper
triangle for a single molecule, reuses coordinates and form factors from
shared memory, and transforms the resulting histogram to reciprocal space on
the GPU. CUDA buffers are retained between calculations. It requires the CUDA
toolkit and an NVIDIA GPU; the CPU implementation remains the default.

```sh
make clean
make GPU=1 CUDA_ARCH=sm_80
./foxs --gpu structure.pdb profile.dat
```

Run the CPU regression and CUDA parity test on a CUDA host with:

```sh
make GPU=1 CUDA_ARCH=sm_80 test
```

`CUDA_ARCH` is optional. Choose the architecture supported by the target GPU.
The `--gpu` flag reports an error when the executable was built without CUDA or
when no CUDA device is available.
