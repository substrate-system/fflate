# Dist Build Restructure Implementation Plan

## Phase 5: Node test suite migration

**Goal:** Repoint the node benchmark harness at `dist/node/index.js` and
move the three empty stub suites out of the node runner.

**Architecture:** The harness keeps spawning CommonJS workers with
`{ eval: true }`. Comparison libraries continue to load with `require`.
Only fflate changes: it is now an ES module, so the worker loads it with
a dynamic `import()` of an absolute `file:` URL.

**Tech Stack:** tsx, @substrate-system/tapzero, node:worker_threads.

**Scope:** Phase 5 of 7.

**Codebase verified:** 2026-08-06

**Depends on:** Phases 3 and 4. The harness loads `dist/node/index.js`,
so a build must have run.

---

## Acceptance Criteria Coverage

This phase implements and tests:

### dist-build-restructure-design.AC5: The node suite passes against the new build

- **dist-build-restructure-design.AC5.1 Success:** `npm run test:node`
  exits 0 with every assertion passing, loading fflate from
  `dist/node/index.js`.
- **dist-build-restructure-design.AC5.2 Success:** The benchmark workers
  spawn and return results for the six fflate methods exercised by the
  node harness (`deflate`, `inflate`, `gzip`, `gunzip`, `zlib`, `unzlib`)
  plus the pako, uzip, tiny-inflate, and zlib comparisons.

  `zip` and `unzip` are not exercised here. The skip at
  `test/2-perf.ts:26` is pre-existing upstream and unchanged by this
  branch, so this phase does not need to justify it -- it only needs to
  avoid claiming coverage it does not have. `unzip` additionally cannot
  round trip through `wc()`: `unzipSync` returns an `Unzipped` plain
  object with no `.buffer`, so the harness's fixed
  `postMessage(buf, [buf.buffer])` transfer list rejects it with
  `DataCloneError`. ZIP coverage lands in Phase 6's browser suite.
- **dist-build-restructure-design.AC5.3 Success:** `test/index.ts`
  imports no ZIP, stream, or one-shot async suite. Those moved to
  `test/browser/`, where they run against the real browser build
  including the minified bundle. Reworded during the final review: the
  original "imports only `./0-valid.js`, `./1-size.js`, and
  `./2-perf.js`" was falsified when `test/3-node-min.ts` was added to
  cover `dist/node/index.min.js`, which the browser suite cannot reach.
  See `test-requirements.md` AC5.3 for the full rationale. The
  three-import code block in Task 2 below is left as sequenced history
  and is correct as of when that task ran.

---

## Investigation findings

Verified against the working tree on 2026-08-06, plus three worker
strategies measured directly against Node 25.8.2.

Confirmed as the design describes:

- `test/util.ts:167` is
  `const fflate = resolve(here, '..', 'lib', 'index.cjs');`, with the
  comment above it explaining that workers `require` it by absolute path
  to bypass the `exports` map.
- `cws()` at `test/util.ts:110-122` builds worker source that calls
  `require()` for both the package under test and `worker_threads`.
- `wc()` at `test/util.ts:129-163` spawns with
  `new Worker(str, { eval: true, workerData, transferList })`.
- `test/index.ts` imports all six suites, `./0-valid.js` through
  `./5-async.js`.
- `test/3-zip.ts`, `test/4-streams.ts`, and `test/5-async.ts` are
  single-line TODO comments with no test content.
- The runner is `@substrate-system/tapzero`, whose API is
  `test(name, async t => { ... })` with `t.ok(value, message)`.

Corrections the design got wrong:

1. **The design's replacement for `wc()` does not work.** The design
   states the generator should emit "ES module source and spawn workers
   via a `data:text/javascript` URL rather than `{ eval: true }`". This
   was measured against Node 25.8.2 and fails. A data: URL module worker
   cannot resolve bare specifiers:

   | strategy | result |
   | --- | --- |
   | `data:` URL module worker, `import pako from 'pako'` | `ERR_UNSUPPORTED_RESOLVE_REQUEST` |
   | `{ eval: true }` CommonJS worker, `require('pako')` | works |
   | `{ eval: true }` CommonJS worker, `import()` of an absolute `file:` URL | works |

   Node's ESM resolver does not resolve bare specifiers relative to a
   `data:` URL, because `data:` is not a special scheme. Every comparison
   library is loaded by bare specifier, so the design's approach breaks
   the entire benchmark harness rather than just the fflate entry.

   Task 1 therefore keeps `{ eval: true }` and changes only how fflate is
   loaded. This is both smaller than the design's proposal and the only
   one of the two that works. `{ eval: true }` remains CommonJS-only in
   current Node, which is fine: the worker body is generated source, not
   a published artifact, and nothing in it needs to be an ES module.

2. **The design's premise for the rework was that "`worker_threads` with
   `{ eval: true }` evaluates as CommonJS, so dropping CJS breaks the
   harness at its foundation."** Only half of that holds. Dropping the
   CJS *build* breaks the harness, because the worker can no longer
   `require` fflate. But `{ eval: true }` itself is not the problem: a
   CommonJS worker can load an ES module through dynamic `import()`,
   which is what Task 1 does. The foundation is intact.

3. **`tiny-inflate` is loaded through the `_cjsDefault` path**
   (`test/util.ts:110`, `wc('tiny-inflate')` with no method), meaning the
   module's export *is* the function. The rewritten accessor must keep
   working for that shape as well as for named exports.

Additional finding:

4. **The design lists jszip as a comparison library, but the harness has
   no jszip worker.** `test/util.ts:169-203` defines workers for fflate,
   pako, uzip, tinyInflate, and zlib only. `jszip` is a devDependency but
   is not wired into `workers`. This phase does not add one: the design's
   Tests section describes preserving existing coverage, not extending
   it. Task 3's expected library list reflects the harness as it is.

5. **`@substrate-system/tapzero` is node-oriented and declares no browser
   export condition.** Its `exports` map has only `import` and `require`
   keys. This does not affect this phase, but it is recorded because
   Phase 6 bundles the browser suite with esbuild, which will resolve the
   `import` condition and inline it. That works, but it is worth knowing
   the package was not explicitly built for a browser target.

---

## Carried forward from the Phase 4 review

`test/util.ts:167` builds the CommonJS path from separate segments:

```ts
const fflate = resolve(here, '..', 'lib', 'index.cjs');
```

Two consequences established on 2026-08-06 and worth having in hand
before this phase starts.

- A `grep` for `lib/` does not find it, because `lib` is its own path
  segment. Search for `'lib'` or for `index.cjs`.
- `npm run test:node` is currently green only because the untracked
  `lib/` directory left over from the old pipeline still exists. Moving
  `lib/` aside makes the suite abort outright with
  `Error: Cannot find module '.../lib/index.cjs' at [worker eval]:2:29`,
  rather than degrade. It would therefore fail on a clean checkout, in
  CI, and as soon as Phase 7 deletes `lib/`.

The replacement target must be loadable from a `worker_threads` eval
string. That does NOT mean it has to be CommonJS. An eval worker is
evaluated as CJS, so bare `require()` of `dist/node/index.js` would
indeed fail -- but a CJS worker can reach an ES module through dynamic
`import()` of an absolute `file:` URL. That is what Task 1 does, and it
was measured working against `dist/node/index.js` on 2026-08-06:
a `{ eval: true }` worker importing the bundle by `file:` URL and
calling `gzipSync` returns correct output.

Verify any fix with `lib/` moved aside.

---

<!-- START_SUBCOMPONENT_A (tasks 1-2) -->

<!-- START_TASK_1 -->
### Task 1: Load fflate from the ESM node bundle in workers

**Verifies:** dist-build-restructure-design.AC5.1,
dist-build-restructure-design.AC5.2

**Files:**
- Modify: `test/util.ts` (the `cws` generator)
- Modify: `test/util.ts` (the `Workerized` declarations through the end
  of the `wc` worker creator)
- Modify: `test/util.ts` (the fflate comment and constant through the
  end of the `workers` object)
- Modify: `test/0-valid.ts` (add the error-marshalling test)

These are deliberately named rather than numbered. Line numbers in this
task would have to be pre-change numbers, since that is the state your
file is in when you start, but every earlier draft of this task drifted
into post-change numbers and one of them pointed at the middle of a
different function. Anchor on the declarations, not on digits.

**Implementation:**

Two coordinated changes. The fflate path moves from the deleted
`lib/index.cjs` to the ESM bundle, expressed as a `file:` URL because it
will be handed to dynamic `import()`. The worker source generator learns
to load `file:` URLs with `import()` while continuing to `require` bare
specifiers.

Replace the `cws` function in `test/util.ts` -- its `// create worker
string` comment through its closing brace -- with:

```ts
// create worker string
//
// Workers stay CommonJS: they are spawned with { eval: true }, which
// Node evaluates as CJS. That is deliberate. A data: URL module worker
// cannot resolve bare specifiers (ERR_UNSUPPORTED_RESOLVE_REQUEST), and
// every comparison library below is loaded by bare specifier.
//
// fflate is now an ES module, so it is loaded with dynamic import() of
// an absolute file: URL. A CommonJS worker can do that.
const cws = (pkg:string, method:string = '_default') => {
  const load = pkg.startsWith('file:') ?
      `await import(${JSON.stringify(pkg)})` :
      `require(${JSON.stringify(pkg)})`;

  // Normalises three shapes: an ESM namespace with named exports, a
  // CJS module whose export is the function, and an ESM default.
  const target = method === '_default' ?
      '(_m.default ?? _m)' :
      `(_m.default ?? _m).${method}`;

  return `
    const { parentPort, workerData } = require('worker_threads');
    (async () => {
      try {
        const _m = ${load};
        const args = Array.isArray(workerData) ?
            workerData :
            [workerData];
        const buf = ${target}(...args);
        parentPort.postMessage(buf, [buf.buffer]);
      } catch (err) {
        const errPayload = err instanceof Error ?
            { name: err.name, message: err.message, stack: err.stack } :
            { name: 'Error', message: String(err) };
        parentPort.postMessage({ err: errPayload });
      }
    })();
  `;
}
```

Then replace the block running from the `export type Workerized`
declaration through the closing brace of `wc` with:

```ts
export type Workerized = (workerData: Uint8Array | [Uint8Array, {}], transferable?: ArrayBuffer[]) => WorkerizedResult;
export interface WorkerizedResult extends PromiseLike<Uint8Array<ArrayBuffer>> {
  timeout(ms: number): void;
};

// Worker creator
const wc = (pkg: string, method?: string): Workerized => {
  const str = cws(pkg, method);
  return (workerData, transferable) => {
    const worker = new Worker(str, {
      eval: true,
      workerData,
      transferList: transferable
    });
    let terminated = false;
    return {
      timeout(ms: number) {
        const tm = setTimeout(() => {
          worker.terminate();
          terminated = true;
        }, ms);
        worker.once('message', () => clearTimeout(tm));
      },
      then(res, rej) {
        return new Promise((res, rej) => {
          worker
            .once('message', msg => {
              if (msg.err) {
                return rej(Object.assign(
                    new Error(msg.err.message),
                    msg.err
                ));
              }
              res(msg);
            })
            .once('error', rej)
            .once('exit', code => {
              if (terminated) rej(new Error('Timed out'));
              else if (code !== 0) rej(new Error('Exited with status code ' + code));
            });
        }).then(res, rej);
      }
    };
  }
}
```

The message handler constructs an Error instance from the structured payload
sent by the worker. It assigns the preserved `{ name, message, stack }` fields
onto the Error object, ensuring callers receive a proper Error instance with
the worker's original stack trace.

Then replace the comment above the fflate constant, the constant itself,
and the whole `workers` object through its closing brace with:

```ts
// Workers load this by absolute URL, which bypasses the "exports" map.
// It is an ES module now, so cws() reaches it with dynamic import().
const fflate = pathToFileURL(
  resolve(here, '..', 'dist', 'node', 'index.js')
).href;

export const workers = {
  fflate: {
    deflate: wc(fflate, 'deflateSync'),
    inflate: wc(fflate, 'inflateSync'),
    gzip: wc(fflate, 'gzipSync'),
    gunzip: wc(fflate, 'gunzipSync'),
    zlib: wc(fflate, 'zlibSync'),
    unzlib: wc(fflate, 'unzlibSync'),
    // zip and unzip are not exercised by the node harness. The skip at
    // test/2-perf.ts:26 is pre-existing upstream and unchanged by this
    // branch.
    //
    // unzip additionally cannot round trip through wc(): it returns an
    // Unzipped plain object with no .buffer, so this harness's fixed
    // postMessage(buf, [buf.buffer]) transfer list rejects it with
    // DataCloneError.
    //
    // ZIP coverage lands in Phase 6's browser suite.
    zip: wc(fflate, 'zipSync'),
    unzip: wc(fflate, 'unzipSync')
  },
  pako: {
    deflate: wc('pako', 'deflateRaw'),
    inflate: wc('pako', 'inflateRaw'),
    gzip: wc('pako', 'gzip'),
    gunzip: wc('pako', 'ungzip'),
    zlib: wc('pako', 'deflate'),
    unzlib: wc('pako', 'inflate')
  },
  uzip: {
    deflate: wc('uzip', 'deflateRaw'),
    inflate: wc('uzip', 'inflateRaw')
  },
  tinyInflate: {
    inflate: wc('tiny-inflate')
  },
  zlib: {
    deflate: wc('zlib', 'deflateRawSync'),
    inflate: wc('zlib', 'inflateRawSync'),
    gzip: wc('zlib', 'gzipSync'),
    gunzip: wc('zlib', 'gunzipSync'),
    zlib: wc('zlib', 'deflateSync'),
    unzlib: wc('zlib', 'inflateSync')
  }
};
```

Add `pathToFileURL` to the existing `url` import on `test/util.ts:3`:

```ts
import { fileURLToPath, pathToFileURL } from 'url';
```

Finally, cover the error-marshalling path. Nothing else in the node
suite reaches either branch of it, and uncovered harness code is how
three review cycles' worth of wrong claims about this file survived a
green run. Add to `test/0-valid.ts`, importing `test` from
`@substrate-system/tapzero` alongside the existing `./util.js` import:

```ts
// Test error-marshalling path: worker throws, structured error is
// reconstructed as Error instance with stack trace
test('worker error handling', async t => {
  const badData = new Uint8Array([9, 9, 9, 9, 9, 9, 9, 9]);
  const promise = workers.fflate.inflate(badData, [badData.buffer]);
  promise.timeout(5000);
  try {
    await promise;
    t.ok(false, 'should have thrown');
  } catch (e) {
    t.ok(e instanceof Error, 'error should be Error instance');
    t.ok(
      typeof e.message === 'string' && e.message.length > 0,
      'error should have non-empty message string'
    );
    // Not just "is a string" -- that is true of any Error. The point of
    // the marshalling is that the stack comes from the worker, so it
    // must not have been constructed locally in wc().
    t.ok(
      typeof e.stack === 'string' && !e.stack.includes('test/util.ts'),
      'error stack should come from the worker, not from wc()'
    );
  }
});
```

Assert on error SHAPE, never on the exact message text: `unexpected EOF`
comes from fflate's internals and is not ours to pin. The third
assertion was verified to bite by mutation -- reducing the handler to
`rej(new Error(msg.err.message))` turns it into `not ok 13`.

**Verification:**

The build must exist. Run `npm run build` first.

Run: `grep -n "lib', 'index.cjs\|lib/index.cjs" test/util.ts`
Expected: no output.

Run: `grep -n "pathToFileURL" test/util.ts`
Expected: two matches, the import and the fflate constant.

Run: `npx tsc --noEmit --project tsconfig.json 2>&1 | grep "error TS" | grep -c "test/util.ts" || true`
Expected: `0`.

Run: `npx tsx test/0-valid.ts`
Expected: exits 0. Every assertion passes. This exercises the fflate
worker (dynamic import of the ESM bundle) against the zlib worker
(`require` of a builtin), so it proves both loader paths at once.

**Commit:** `test: load fflate esm bundle in benchmark workers`
<!-- END_TASK_1 -->

<!-- START_TASK_2 -->
### Task 2: Move the stub suites out of the node runner

**Verifies:** dist-build-restructure-design.AC5.3

**Files:**
- Move: `test/3-zip.ts` to `test/browser/zip.ts`
- Move: `test/4-streams.ts` to `test/browser/streams.ts`
- Move: `test/5-async.ts` to `test/browser/async.ts`
- Modify: `test/index.ts`

**Implementation:**

The three files are single-line TODO comments. They move rather than
being deleted and recreated, so the intent recorded in each comment
follows into the browser suite where Phase 6 implements it.

```bash
mkdir -p test/browser
git mv test/3-zip.ts test/browser/zip.ts
git mv test/4-streams.ts test/browser/streams.ts
git mv test/5-async.ts test/browser/async.ts
```

Then replace `test/index.ts` with:

```ts
// Entry point for the node test suite.
//
// tapzero has no CLI runner: importing a suite registers its tests, and
// the runner flushes them in registration order once this module
// finishes evaluating. The numeric filename prefixes therefore set the
// run order.
//
// ZIP, stream, and async coverage lives in test/browser/, where it runs
// against the real browser build including the minified bundle.
import './0-valid.js';
import './1-size.js';
import './2-perf.js';
```

**Verification:**

Run: `ls test`
Expected: `0-valid.ts`, `1-size.ts`, `2-perf.ts`, `browser`, `data`,
`index.ts`, `results`, `util.ts`. No numbered files 3 through 5.

Run: `ls test/browser`
Expected: `async.ts`, `streams.ts`, `zip.ts`.

Run: `grep -cE "3-zip|4-streams|5-async" test/index.ts || true`
Expected: `0`.

**Commit:** `test: move stub suites to browser directory`
<!-- END_TASK_2 -->

<!-- END_SUBCOMPONENT_A -->

<!-- START_TASK_3 -->
### Task 3: Verify the full node suite

**Verifies:** dist-build-restructure-design.AC5.1,
dist-build-restructure-design.AC5.2,
dist-build-restructure-design.AC5.3

**Files:**
- No files changed. This task is a gate.

**Implementation:**

Run the whole node suite against a fresh build.

Be aware that `test/util.ts:13-24` downloads fixtures from the public
internet on first run, caching them under `test/data/`. The first
execution is slow and requires network access. `test/1-size.ts` and
`test/2-perf.ts` are benchmarks over multi-megabyte inputs and take
noticeably longer than `test/0-valid.ts`.

**Verification:**

Run: `npm run build && npm run test:node`
Expected: exits 0. TAP output ends with a plan line and no `not ok`
entries.

Run: `npm run test:node 2>&1 | grep -c "^not ok" || true`
Expected: `0`.

Run: `npm run test:node 2>&1 | grep -cE "ERR_UNSUPPORTED_RESOLVE_REQUEST|Cannot find module|ERR_MODULE_NOT_FOUND" || true`
Expected: `0`. Any match means a worker failed to load its module.

Run: `test -f test/results/timings.json && test -f test/results/longTimings.json && echo "results written"`
Expected: prints `results written`. This confirms the benchmark suites
ran to completion rather than being skipped.

Run: `node -e "const t=require('./test/results/timings.json'); const libs=new Set(Object.keys(t).map(k=>k.split('.')[0])); console.log([...libs].sort().join(','))"`
Expected: `fflate,pako,tinyInflate,uzip,zlib`. A missing library means
its worker never produced a result.

**Commit:** `test: verify node suite against dist build`
<!-- END_TASK_3 -->

---

## Phase 5 completion criteria

- `npm run test:node` exits 0 with zero `not ok` lines.
- `test/util.ts` references `dist/node/index.js` as a `file:` URL and no
  longer mentions `lib/index.cjs`.
- Benchmark results are written for all five libraries.
- `test/index.ts` imports exactly three suites.
- `test/browser/` holds the three moved stubs.

`npm test` still fails at `test:browser`, which Phase 6 implements.
