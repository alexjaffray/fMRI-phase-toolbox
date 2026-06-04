# History and Roadmap

## Legacy v1 Content

The original project contained the following:

- hard-coded local paths for data, QSM, MRecon, and outputs

With spaghetti code and scripts performing:

- interactive file selection
- image loading and rotation
- QSM phase unwrapping and background-field removal
- SVD component selection
- solid-harmonic regression
- physiologic log comparison
- plotting and manuscript output generation

## Goals in v2

1. Converge on a MATLAB API and interface
2. Keep input/output explicit through a configuration struct
3. Keep algorithmic units small enough to test 
5. Separate third-party/vendor code from original toolbox code 
