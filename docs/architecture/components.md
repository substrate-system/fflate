# Components

Source modules that ship, and how the async APIs get work onto a worker.

Sources:

- [`src/index.ts`](../../src/index.ts) -- every public export, plus the
  compression internals and the worker plumbing.
- [`src/worker.ts`](../../src/worker.ts) -- browser worker spawner.
- [`src/node-worker.ts`](../../src/node-worker.ts) -- node worker spawner.
- [`scripts/build.ts`](../../scripts/build.ts) -- performs the swap between
  the two spawners.

## Modules

| Module | Ships as | Role |
| --- | --- | --- |
| `src/index.ts` | Bundled into every artifact | Sole public entry point. Sync codecs, streaming classes, async wrappers, ZIP/gzip/zlib containers. ~4200 lines. |
| `src/node-worker.ts` | Bundled into `dist/node/` | Spawns a `node:worker_threads` `Worker` with `eval:true`. |
| `src/worker.ts` | Bundled into `dist/browser/` and `dist/umd/` | Spawns a DOM `Worker` from a `Blob` object URL. |

`src/index.ts` imports `./node-worker` unconditionally. The browser and UMD
builds rewrite that specifier to `./worker.ts` through an esbuild `onResolve`
plugin named `browser-worker-swap`. `--alias` cannot do this, because it only
rewrites bare specifiers and this one is relative.

Neither worker module is a separate published entry point. Only
`index.d.ts` is emitted for consumers; the build deletes `worker.d.ts` and
`node-worker.d.ts` after declaration emit.

## Async execution path

The async APIs do not load a worker script from a URL. They reconstruct the
worker's source as a string at call time, from the functions already in the
bundle, and hand that string to the platform spawner.

The chain, all in `src/index.ts`:

1. `bInflt` / `bDflt` (lines 1153-1154) return arrays of the primitives an
   inflate or deflate worker needs. Container-specific extras come from
   `gze`, `guze`, `zle`, `zule` (lines 1157-1163).
2. `wcln` (line 1098) serializes one such array. It calls the function to get
   the values, then calls `.toString()` on the function and slices the text
   between the first `[` and the last `]` to recover the identifier **names**.
   Values pair with names positionally.
3. `wrkr` (line 1142) concatenates the serialized groups, appends an
   `onmessage` bootstrap, and calls the spawner. It also wraps the caller's
   callback in `rhyd`, the single point where a worker error regains its
   `FlateError` prototype (see "Spawner contracts" below).
4. `cbify` (line 1175) wraps one-shot calls; `astrmify` (line 1214) wraps the
   streaming classes.

Step 2 is the load-bearing oddity. Name recovery reads source text, so
anything that rewrites the text of those array literals -- a minifier
inlining a constant, or a bundler wrapping functions -- breaks it. The failure
is silent: sync APIs never touch this path and keep passing. See
[build-and-packaging.md](build-and-packaging.md) for the flags this
constrains.

`wcln` handles a class the same way it handles a function, emitting a named
class expression, which is why `FlateError` can sit in `bInflt` alongside
`err` (which constructs it inside the worker). One caveat: the branch that
copies prototype members uses `for (const t in v.prototype)`, and class
methods are non-enumerable. A class with methods would have them silently
dropped from the worker copy. `FlateError` has only a constructor, which
`toString()` already includes.

`bDflt` carries no `err` and no `FlateError`: no deflate-side code path
constructs one.

## Spawner contracts

Both spawners share the signature
`(code, id, msg, transfer, cb) => Worker`.

**Browser** (`src/worker.ts`) appends an `error` listener that forwards
`[message, code, stack]` under a `$e$` key, then rebuilds an `Error` on the
main thread. It caches one object URL per `id` in a module-level `ch2` map.

**Node** (`src/node-worker.ts`) appends a prelude that bridges
`parentPort` to the `onmessage` / `postMessage` globals the shared code
expects, and aliases `self` to `global`. It filters untransferable buffers
out of the transfer list via `isMarkedAsUntransferable`, and overrides
`terminate` so a deliberate termination does not report as a nonzero exit.

Neither spawner can deliver a `FlateError` as such. The browser rebuilds a
plain `Error` from `[message, code, stack]`; node re-serializes the error
for its `error` event. Both preserve `message`, `stack`, `name`, and the
numeric `code`, but neither preserves the prototype. `rhyd` in
`src/index.ts` reconstructs a real `FlateError` from any error carrying a
numeric `code`, keeping the worker-side stack. Errors without one -- a
worker exiting nonzero, an out-of-memory -- pass through untouched, so an
unrelated failure is never mislabelled as a corrupt-archive error. Because
this sits in `wrkr`, every async API and every worker-backed streaming
class gets it; do not add error handling that bypasses `wrkr`.

## Worker cache ids

Each async API passes a fixed integer id, which keys the serialized-source
cache described in [runtime-state.md](runtime-state.md).

| Id | API | Kind |
| --- | --- | --- |
| 0 | `deflate` | one-shot |
| 1 | `inflate` | one-shot |
| 2 | `gzip` | one-shot |
| 3 | `gunzip` | one-shot |
| 4 | `zlib` | one-shot |
| 5 | `unzlib` | one-shot |
| 6 | `AsyncDeflate` | streaming |
| 7 | `AsyncInflate` | streaming |
| 8 | `AsyncGzip` | streaming |
| 9 | `AsyncGunzip` | streaming |
| 10 | `AsyncZlib` | streaming |
| 11 | `AsyncUnzlib` | streaming |

`AsyncUnzipInflate` (line 3910) has no id of its own. Below 320000 bytes it
uses the synchronous `Inflate`; at or above that it delegates to
`AsyncInflate`, so it shares id 7.
