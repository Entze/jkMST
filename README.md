# jkMST

jkMST is a Julia Program to optimize k-MST problems.

## Dependencies

- Julia v1.4+
- CPLEX 12.10 or 21.1

## Installation

To install the script locally you can run:

```bash
make
```
The script attempts to find the CPLEX installation in a set of common locations, if it does not find the installation set the environment variable `CPLEX_STUDIO_BINARIES`:

```bash
export CPLEX_STUDIO_BINARIES=/path/to/cplex
```

## Usage

After installing the bash script can be run.

For example run

```bash
./jkMST --help
```

for details.

The most common command to run will be:

```bash
./jkMST --directory data # finds all .dat files in all subdirectories
    --mode mtz scf mcf # launches every instance in every mode
    --size 0.2 0.5 # relative sizes of the k-MST tree
    --timeout 600 # 10 minute runtime in seconds
    --generate-timeout 60 # 1 minute time for building the model
```

The flag `--verbose` outputs quite a bit of additional information which might be interesting.

## Contributing

This is university project and therefore it is highly academic/non-practical in nature.

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)