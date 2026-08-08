# Runtime state

What this library holds in memory across calls, and what it deliberately
does not hold at all.

Sources:

- [`src/index.ts`](../../src/index.ts) -- `ch` cache, `wrkr`, `cbfs`.
- [`src/worker.ts`](../../src/worker.ts) -- `ch2` object URL cache.

## No durable state

There is no persisted state of any kind. The package writes no files, opens
no connections, reads no environment variables, and keeps nothing across
process boundaries. Nothing here has a TTL, a backup story, or a security
boundary beyond the calling process. A consumer that unloads the module
loses everything below and pays only a recomputation cost.

## Process-local caches

Both caches are module-level, unbounded in principle but bounded in practice
by the fixed id set in [components.md](components.md). Neither is ever
evicted or invalidated.

| Name | Location | Key | Value | Lifetime |
| --- | --- | --- | --- | --- |
| `ch` | `src/index.ts:1129` | worker id, `0`-`11` | `{ c:string, e:Record<string, unknown> }` -- serialized worker source and the non-function values pulled out of it | Module lifetime |
| `ch2` | `src/worker.ts:1` | worker id, as a string key | Blob object URL | Module lifetime, browser builds only |

`ch` is populated lazily on the first async call for a given id, inside
`wrkr` (`src/index.ts:1143`). Serializing the source is the expensive part,
so every later call for that id reuses the string.

`ch2` is browser-only. `URL.createObjectURL` is called once per id and the
URL is never revoked, so the underlying Blob stays alive for the life of the
document. That is intentional -- revoking it would defeat the cache -- but it
does mean object URLs are a small, fixed, permanent allocation once the async
APIs are used.

The node spawner has no equivalent cache. It passes the source string
straight to `new Worker(..., { eval:true })` on every call.

## Per-call state

`ch[id].e` holds the typed-array constants a worker needs. These cannot be
shared across workers, because posting them transfers ownership of their
buffers. `wrkr` therefore copies the record with `mrg` and passes it through
`cbfs` (`src/index.ts:1131`), which reallocates each typed array and returns
the fresh buffers as the transfer list. The cached originals are never
transferred and stay intact for the next call.

Worker instances themselves are not pooled. Each async call or streaming
instance owns its worker and is responsible for terminating it, which is
what the `terminate` member on the async classes and the `AsyncTerminable`
return values are for.
