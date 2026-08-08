# Architecture

The architecture inventory lives in [`architecture/`](architecture/). Start
at [`architecture/INDEX.md`](architecture/INDEX.md).

- [Components](architecture/components.md) -- source modules, async worker
  execution, the browser/node swap.
- [Runtime state](architecture/runtime-state.md) -- process-local caches;
  nothing is persisted.
- [Build and packaging](architecture/build-and-packaging.md) -- the pipeline,
  the `dist/` artifacts, the `exports` map, and the silent invariants.

For the public API, see the generated reference in [`README.md`](README.md)
and the `classes/`, `functions/`, and `interfaces/` directories. For design
rationale, see [`../specs/`](../specs/).
