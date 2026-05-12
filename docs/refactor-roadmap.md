# Refactor Roadmap

## What v1 Contains

The original project combines several concerns in research scripts:

- hard-coded local paths for data, QSM, MRecon, and outputs
- interactive file selection
- image loading and rotation
- QSM phase unwrapping and background-field removal
- SVD component selection
- solid-harmonic regression
- physiologic log comparison
- plotting and manuscript output generation

## v2 Direction

1. Stabilize a MATLAB API first.
2. Keep input/output explicit through a configuration struct.
3. Keep algorithmic units small enough to test with synthetic arrays.
4. Preserve original scripts under `legacy/` while refactoring behavior.
5. Separate third-party/vendor code from original toolbox code.
6. Add Python or Julia only after MATLAB behavior is validated.

## Candidate Language Strategy

MATLAB should remain the reference implementation because the current workflow
depends on MATLAB-native NIfTI, Image Processing Toolbox operations, QSM.m, and
MRecon. Python is the likely best second target for community adoption because
the neuroimaging ecosystem already has NiBabel, Nilearn, NumPy, SciPy, and
BIDS tooling. Julia may be useful for performance-oriented numerical kernels,
but should probably follow once the MATLAB API and test fixtures are stable.

## Near-Term Tasks

- Add synthetic tests for harmonic basis projection.
- Add synthetic tests for SVD component selection.
- Add a small public example dataset or documented fixture generation.
- Replace hard-coded protocol assumptions such as phase scaling with config.
- Define BIDS-friendly input discovery.
- Add continuous integration for MATLAB tests if the repository is published.
