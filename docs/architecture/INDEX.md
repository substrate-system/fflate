# Architecture Inventory

Current runtime facts for `@substrate-system/fflate`: what exists, where it
lives, and the contract it holds today. Rationale and alternatives belong in
design documents under `specs/`; this tree records the runtime as built.

Last verified: 2026-08-08, against commit `b5ed05b` on branch `issues`.

## Categories

| File | Covers |
| --- | --- |
| [components.md](components.md) | Source modules, the async worker execution paths, and the platform swap that picks between them. |
| [runtime-state.md](runtime-state.md) | Process-local state held by the library at runtime, and what is deliberately absent. |
| [build-and-packaging.md](build-and-packaging.md) | The build pipeline, the artifacts in `dist/`, the `exports` map, and the invariants that are silent when broken. |

## Scope notes

This package is a dependency-free ESM library. It opens no sockets, reads no
configuration, and persists nothing between processes. Categories that a
network service would carry -- durable stores, subjects, endpoints, auth
policies, deployment topology -- have no counterpart here and are not stubbed
out.

There is no ADR or FDR tree in this repository, and no `docs/GLOSSARY.md`.
Design rationale currently lives in `specs/` and in per-change plans under
`docs/implementation-plans/`. Inventory files link to those instead.

## Relationship to generated docs

`docs/` also holds typedoc output (`docs/classes/`, `docs/functions/`,
`docs/interfaces/`, `docs/type-aliases/`, `docs/variables/`, and
`docs/README.md`). That output is generated and must not be hand-edited. This
`architecture/` subtree is hand-maintained and survives regeneration because
`build-docs` runs with `--cleanOutputDir false`
(`package.json`, `scripts.build-docs`). Keep that flag.
