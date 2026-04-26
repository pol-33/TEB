# TEB Sequence Mapper Project

This repository contains our two implementations for the TEB sequence-mapper project:

- [mapper-speed](mapper-speed/README.md): the speed-focused version, optimized for throughput.
- [mapper-memory](mapper-memory/README.md): the memory-focused version, optimized for lower index size and lower peak RSS.

## Authors

- Pol Plana Torrents
- Joan Teruel
- Mariam Delgado

## Project Structure

- [mapper-speed](mapper-speed/README.md): short usage guide for the speed-oriented mapper.
- [mapper-speed documentation](mapper-speed/DOCUMENTATION.md): full design and benchmarking notes for the speed-oriented mapper.
- [mapper-memory](mapper-memory/README.md): short usage guide for the memory-oriented mapper.
- [mapper-memory documentation](mapper-memory/DOCUMENTATION.md): full design and benchmarking notes for the memory-oriented mapper.
- [data](data/README.md): notes about dataset placement and local inputs.
- `legacy/`: preserved earlier coursework code.

## Build

Build the speed-oriented mapper:

```bash
cd mapper-speed
make
```

Build the memory-oriented mapper:

```bash
cd mapper-memory
make
```

## Choosing a Version

Use `mapper-speed` if:

- you want the best runtime,
- you are benchmarking for throughput,
- you want the implementation intended for the speed track.

Use `mapper-memory` if:

- you want the smallest / lowest-RSS index configuration,
- you are benchmarking for the memory track,
- you are willing to trade a lot of runtime for lower memory usage.
