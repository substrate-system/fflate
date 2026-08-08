# Unzip mtime Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps
> use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose each ZIP entry's last-modification time as a decoded
`Date` on `UnzipFileInfo` (the `unzip`/`unzipSync` filter) and on
`UnzipFile` (the streaming `Unzip` class).

**Architecture:** ZIP headers store modification time as one 32-bit
little-endian DOS date-time word -- at offset `+12` in a central
directory header, `+10` in a local file header. A single internal helper,
`dosdt`, unpacks that word into a `Date`. `zh()` (the central-header
parser) returns the raw word as a new tuple element; the streaming
`Unzip` class reads the local header itself. Both then call `dosdt`.

**Tech Stack:** TypeScript (`target: es2022`), esbuild + terser, tsx,
tapzero assertions, tapout browser runner.

**Spec:** [`specs/2026-08-08-unzip-mtime-design.md`](../../specs/2026-08-08-unzip-mtime-design.md)

## Global Constraints

- **One npm invocation at a time.** Every script shares `dist/`. A
  concurrent run produces an `npm test` exit 1 that looks like a real
  failure and is not.
- **`npm run test:browser` requires a prior `npm run build`.**
  `test/browser/index.ts` imports from `../../dist/browser/index.js`,
  which is gitignored. Always run `npm run build` first.
- **Never run** `npm run gh-pages` or `npm run build-docs`. The first is
  destructive to the `gh-pages` branch; the second overwrites `docs/`,
  which also holds this plan.
- **Do not hand-edit `docs/interfaces/`, `docs/classes/`, or any other
  typedoc output.** It regenerates. `docs/implementation-plans/` and
  `docs/test-plans/` are hand-written and are the exception.
- **Never verify anything with `node --input-type=module -e`.** That form
  leaks `--input-type=module` into `process.execArgv`, worker threads
  inherit it, and fflate's worker payload then throws
  `ReferenceError: u8 is not defined`. Use a real `.mjs`/`.cjs` file.
- **Max line length 80 columns.** No space between colon and type
  annotation in TypeScript (`mtime:Date`, not `mtime: Date`) **in test
  files**. `src/index.ts` is upstream fflate code and uses the opposite
  convention (`mtime: Date`) -- match the file you are editing.
- **No em dashes and no `->` arrow characters** in Markdown or comments.
  Use `--` and `->` respectively.
- **`keepNames` stays `false`** in every esbuild config, and no new terser
  `compress` options. Neither task here touches `scripts/build.ts`; if you
  find yourself editing it, stop and re-read the spec.

---

## File Structure

| File | Responsibility |
| --- | --- |
| `src/index.ts` | The `dosdt` helper, the `zh` tuple change, the `mtime` members on `UnzipFileInfo`/`UnzipFile`, and the three sites that construct those objects. |
| `test/browser/zip.ts` | All four tests. This is the only file in the project with ZIP coverage, and the browser suite runs it twice -- once against `dist/browser/index.js`, once against `dist/browser/index.min.js`. |
| `README.md` | The two hand-written examples that show unzip metadata. |
| `TODO/README.md` | The issue checkbox. |

`test/0-valid.ts` (the node suite) has no ZIP coverage at all and is not
the right home for these tests. Do not add them there.

---

## Correction, applied during execution

Two claims in the original text of this plan were wrong. Both are fixed
in place below; recorded here so the diff is not mysterious.

1. **Tasks 1 and 2 cannot land separately.** `npm run build` runs
   `tsc --emitDeclarationOnly`, so making `UnzipFileInfo.mtime`
   *required* breaks declaration emit at the async `unzip` filter
   literal until that call site is updated too. Task 1 leaves the tree
   unbuildable on its own. The two were merged into a single commit,
   with Task 2's test written first and the compile error serving as
   its red state.
2. **The 2099 test does not catch a signed shift.** Measured: with the
   `& 127` mask present, `((t >> 25) & 127) + 1980` and
   `((t >>> 25) & 127) + 1980` both return 2099 for `0xEF9FBF7D`. The
   mask discards the sign extension. `>>>` is still correct and kept,
   but for a different reason -- see the decoder comment.

---

## Task 1: The decoder and the `unzipSync` filter

**Files:**
- Modify: `src/index.ts` (helper near :2862, `zh` at :2865,
  `UnzipFileInfo` at :3572, `unzipSync` at :3967)
- Test: `test/browser/zip.ts`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `const dosdt = (t:number) => Date` -- internal, unexported. Tasks 2
    and 3 call it.
  - `zh(d, b, z)` returns a 7-element readonly tuple
    `[compression, compressedSize, originalSize, filename, nextOffset,
    localHeaderOffset, rawDosDateTime]`. Task 2 destructures the 7th.
  - `UnzipFileInfo.mtime:Date` -- non-optional. Task 2 relies on it.
  - Test constants `MTIME_A`, `MTIME_B`, `MTIME_MIN`, `MTIME_MAX` in
    `test/browser/zip.ts`. Tasks 2, 3, and 4 reuse them.

- [ ] **Step 1: Add the shared test constants**

At the top of `test/browser/zip.ts`, directly below the existing
`const MTIME = new Date('2020-01-01T00:00:00Z')` (line 18), add:

```ts
// DOS date-time fixtures. The format stores seconds/2 and carries no
// time zone, so every value here uses an even second, zero
// milliseconds, and the local-time Date constructor -- the same one
// wzh writes with (src/index.ts:2917). Odd seconds or a UTC
// constructor would make these round-trips fail for reasons that have
// nothing to do with the decoder.
const MTIME_A = new Date(2021, 4, 17, 13, 24, 36)
const MTIME_B = new Date(1999, 11, 31, 23, 59, 58)

// The ends of the representable DOS range. wzh throws error 10 for any
// year outside 1980..2099, so these are the true boundaries.
const MTIME_MIN = new Date(1980, 0, 1, 0, 0, 0)
const MTIME_MAX = new Date(2099, 11, 31, 23, 59, 58)
```

- [ ] **Step 2: Write the two failing tests**

Append to the `zipSuite` function body in `test/browser/zip.ts`, before
its closing brace:

```ts
  // The DOS mod-time word lives in the central directory header, which
  // is the header zh() parses for unzipSync. Two entries with
  // DIFFERENT times is the point: a single shared read would pass an
  // equal-times test.
  test(`${label} > unzipSync filter receives per-entry mtime`, async t => {
    const compressible = fx()[0]
    const random = fx()[1]

    const archive = f.zipSync({
      'a.bin': [compressible, { mtime: MTIME_A }],
      'b.bin': [random, { mtime: MTIME_B }]
    })

    const seen:Record<string, number> = {}
    f.unzipSync(archive, {
      filter (file) {
        seen[file.name] = file.mtime.getTime()
        return true
      }
    })

    t.equal(seen['a.bin'], MTIME_A.getTime(),
      'a.bin mtime round-trips through the unzipSync filter')
    t.equal(seen['b.bin'], MTIME_B.getTime(),
      'b.bin mtime round-trips through the unzipSync filter')
  })

  test(`${label} > unzipSync filter decodes DOS mtime boundaries`, async t => {
    const compressible = fx()[0]
    const random = fx()[1]

    const archive = f.zipSync({
      'min.bin': [compressible, { mtime: MTIME_MIN }],
      'max.bin': [random, { mtime: MTIME_MAX }]
    })

    const seen:Record<string, number> = {}
    f.unzipSync(archive, {
      filter (file) {
        seen[file.name] = file.mtime.getTime()
        return true
      }
    })

    t.equal(seen['min.bin'], MTIME_MIN.getTime(),
      'DOS epoch 1980-01-01 round-trips')
    // 2099 packs as 0xEF9FBF7D, with bit 31 set. This pins the top of
    // the DOS range and the year field's shift distance and mask
    // width. It does NOT distinguish `>>` from `>>>` in dosdt: the
    // `& 127` discards a signed shift's sign extension, so both give
    // 2099. Measured, not assumed.
    t.equal(seen['max.bin'], MTIME_MAX.getTime(),
      'DOS maximum 2099-12-31 round-trips')
  })
```

- [ ] **Step 3: Run the tests to verify they fail**

```bash
npm run build && npm run test:browser
```

Expected: FAIL. The four `t.equal` assertions report `not ok`, because
`file.mtime` is `undefined` and `undefined.getTime()` throws, or the
suite reports the `TypeError`. Either way the two new test names must
appear as failures. If they appear as passes, the tests are not wired up
-- stop and fix that before continuing.

- [ ] **Step 4: Add the `dosdt` helper**

In `src/index.ts`, immediately after the `slzh` definition (:2862) and
before `// read zip header`, insert:

```ts
// decode a DOS date-time dword into a Date. The word packs, from the
// top: year-1980 (7 bits), month (4), day (5), hour (5), minute (6),
// seconds/2 (5).
//
// Two things here are easy to get wrong. Shifts are unsigned because
// any year from 2044 on sets bit 31 (2099 packs as 0xEF9FBF7D). The
// `& 127` happens to mask a signed shift's sign extension away too, so
// `>>` would pass the tests as written -- but `(t >>> 25) + 1980` is a
// tempting simplification of that first line, and it is only correct
// with the unsigned shift. And the month mask must be parenthesised
// before the subtraction: `& 15 - 1` binds as `& 14` and silently
// clears the field's low bit.
const dosdt = (t: number) => new Date(
  ((t >>> 25) & 127) + 1980,
  ((t >>> 21) & 15) - 1,
  (t >>> 16) & 31,
  (t >>> 11) & 31,
  (t >>> 5) & 63,
  (t & 31) * 2
);
```

- [ ] **Step 5: Return the raw dword from `zh`**

Replace the return statement of `zh` (`src/index.ts:2868`). Before:

```ts
  return [b2(d, b + 10), sc, su, fn, es + efl + b2(d, b + 32), off] as const;
```

After:

```ts
  // b + 12 is the DOS mod-time dword. It stays RAW here: unzip and
  // unzipSync guard the filter with `!fltr || fltr({...})`, so
  // short-circuiting means no Date is allocated when no filter is
  // supplied. Returning a Date instead would allocate one per central
  // directory entry on every unzip.
  return [
    b2(d, b + 10), sc, su, fn, es + efl + b2(d, b + 32), off, b4(d, b + 12)
  ] as const;
```

- [ ] **Step 6: Add `mtime` to `UnzipFileInfo`**

In `src/index.ts`, inside `interface UnzipFileInfo`, after the
`compression` member (ends :3594), add:

```ts

  /**
   * The last modification time of the file. DOS timestamps carry no
   * time zone and have two-second resolution, so this is a local-time
   * Date on an even second, with no milliseconds.
   */
  mtime: Date;
```

- [ ] **Step 7: Wire it into `unzipSync`**

In `unzipSync` (`src/index.ts:3966-3974`), extend the destructuring and
the filter argument. Before:

```ts
    const [c, sc, su, fn, no, off] = zh(data, o, z), b = slzh(data, off);
    o = no;
    if (!fltr || fltr({
      name: fn,
      size: sc,
      originalSize: su,
      compression: c
    })) {
```

After:

```ts
    const [c, sc, su, fn, no, off, mtm] = zh(data, o, z), b = slzh(data, off);
    o = no;
    if (!fltr || fltr({
      name: fn,
      size: sc,
      originalSize: su,
      compression: c,
      mtime: dosdt(mtm)
    })) {
```

The name is `mtm`, not `mt`. `mt` is already a module-level binding at
`src/index.ts:3852` -- the `queueMicrotask`/`setTimeout` scheduler that
`unzip` calls. Shadowing it inside these loops is a trap waiting for the
next person; do not.

- [ ] **Step 8: Run the tests to verify they pass**

```bash
npm run build && npm run test:browser
```

Expected: PASS. Both new test names appear twice each -- once under the
`plain` label and once under `min` -- and all four assertions are `ok`.

- [ ] **Step 9: Typecheck**

```bash
npx tsc --noEmit -p tsconfig.json > /dev/null; echo "EXIT:$?"
```

Expected: `EXIT:0`. (`tsconfig.json` sets `listFiles: true`, so the
redirect suppresses a long file dump. Only the exit code matters.)

- [ ] **Step 10: Run the full suite**

```bash
npm test
```

Expected: PASS. This also runs the node suite, which must not regress.

- [ ] **Step 11: Commit**

```bash
git add src/index.ts test/browser/zip.ts
git commit -m "feat: decode entry mtime for the unzipSync filter

Adds an internal dosdt helper that unpacks the DOS date-time dword
from a ZIP header into a Date, returns that dword from zh, and
surfaces it as UnzipFileInfo.mtime."
```

---

## Task 2: The async `unzip` filter

**Files:**
- Modify: `src/index.ts` (`unzip` at :3904-3920)
- Test: `test/browser/zip.ts`

**Interfaces:**
- Consumes: `dosdt`, the 7-element `zh` tuple, and
  `UnzipFileInfo.mtime` from Task 1. `MTIME_A` and `MTIME_B` from Task 1.
- Produces: nothing new. `unzip`'s filter now receives the same
  `UnzipFileInfo` shape `unzipSync`'s does.

- [ ] **Step 1: Write the failing test**

Append to the `zipSuite` function body in `test/browser/zip.ts`:

```ts
  // Same central-header path as unzipSync, but through the async
  // control flow. The existing unzipAsync() helper in this file uses
  // the two-argument form, which cannot carry a filter, so this test
  // drives f.unzip directly.
  test(`${label} > unzip filter receives per-entry mtime`, async t => {
    const compressible = fx()[0]
    const random = fx()[1]

    const archive = f.zipSync({
      'a.bin': [compressible, { mtime: MTIME_A }],
      'b.bin': [random, { mtime: MTIME_B }]
    })

    const seen:Record<string, number> = {}
    try {
      await withTimeout(() => new Promise<void>((resolve, reject) => {
        f.unzip(archive, {
          filter (file) {
            seen[file.name] = file.mtime.getTime()
            return true
          }
        }, err => {
          if (err) reject(err)
          else resolve()
        })
      }))

      t.equal(seen['a.bin'], MTIME_A.getTime(),
        'a.bin mtime round-trips through the unzip filter')
      t.equal(seen['b.bin'], MTIME_B.getTime(),
        'b.bin mtime round-trips through the unzip filter')
    } catch (e) {
      t.ok(false, `unzip filter mtime failed: ${e}`)
    }
  })
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
npm run build && npm run test:browser
```

Expected: FAIL. `unzip filter receives per-entry mtime` reports `not ok`
-- `file.mtime` is `undefined` in this code path even though Task 1
landed, because `unzip` builds its own filter object literal.

- [ ] **Step 3: Wire `mtime` into `unzip`**

In `unzip` (`src/index.ts:3904` and :3915-3920). Before:

```ts
      const [c, sc, su, fn, no, off] = zh(data, o, z), b = slzh(data, off);
```

After:

```ts
      const [c, sc, su, fn, no, off, mtm] = zh(data, o, z), b = slzh(data, off);
```

And before:

```ts
      if (!fltr || fltr({
        name: fn,
        size: sc,
        originalSize: su,
        compression: c
      })) {
```

After:

```ts
      if (!fltr || fltr({
        name: fn,
        size: sc,
        originalSize: su,
        compression: c,
        mtime: dosdt(mtm)
      })) {
```

Again: `mtm`, not `mt`. `unzip` itself calls `mt(...)` at :3879 and
:3881 to schedule the callback. Naming the destructured variable `mt`
shadows the scheduler.

- [ ] **Step 4: Run the test to verify it passes**

```bash
npm run build && npm run test:browser
```

Expected: PASS. The new test name appears under both `plain` and `min`,
with both assertions `ok`, and Task 1's tests still pass.

- [ ] **Step 5: Typecheck**

```bash
npx tsc --noEmit -p tsconfig.json > /dev/null; echo "EXIT:$?"
```

Expected: `EXIT:0`.

- [ ] **Step 6: Commit**

```bash
git add src/index.ts test/browser/zip.ts
git commit -m "feat: decode entry mtime for the async unzip filter"
```

---

## Task 3: The streaming `Unzip` class

**Files:**
- Modify: `src/index.ts` (`UnzipFile` at :3607, `Unzip.push` at
  :3776-3808)
- Test: `test/browser/zip.ts`

**Interfaces:**
- Consumes: `dosdt` from Task 1. `MTIME_A` and `MTIME_B` from Task 1.
- Produces: `UnzipFile.mtime:Date` -- non-optional.

This is genuinely separate coverage, not a duplicate of Tasks 1 and 2.
`Unzip` never calls `zh`; it parses the **local** file header, where the
mod-time dword sits at `+10` rather than `+12`.

- [ ] **Step 1: Write the failing test**

Append to the `zipSuite` function body in `test/browser/zip.ts`:

```ts
  // Unzip parses LOCAL file headers, not the central directory, so
  // this exercises a different offset (+10) and a different code path
  // from the two filter tests above.
  test(`${label} > Unzip stream exposes per-entry mtime`, async t => {
    const compressible = fx()[0]
    const random = fx()[1]

    const archive = f.zipSync({
      'a.bin': [compressible, { mtime: MTIME_A }],
      'b.bin': [random, { mtime: MTIME_B }]
    })

    const seen:Record<string, number> = {}
    // No stream is ever started: this test reads header metadata only.
    // The decoder is registered anyway because start() would throw
    // without one, and a later edit that adds a start() call should
    // not have to rediscover that.
    const unzipper = new f.Unzip(file => {
      seen[file.name] = file.mtime.getTime()
    })
    unzipper.register(f.UnzipInflate)
    unzipper.push(archive, true)

    t.equal(seen['a.bin'], MTIME_A.getTime(),
      'a.bin mtime read from its local header')
    t.equal(seen['b.bin'], MTIME_B.getTime(),
      'b.bin mtime read from its local header')
  })
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
npm run build && npm run test:browser
```

Expected: FAIL. `Unzip stream exposes per-entry mtime` reports `not ok`;
`file.mtime` is `undefined` because `Unzip.push` builds its own file
object literal.

- [ ] **Step 3: Add `mtime` to `UnzipFile`**

In `src/index.ts`, inside `interface UnzipFile`, after the
`originalSize?: number;` member (ends :3636), add:

```ts

  /**
   * The last modification time of the file. DOS timestamps carry no
   * time zone and have two-second resolution, so this is a local-time
   * Date on an even second, with no milliseconds. Unlike size and
   * originalSize, this is always present -- the mod-time field is a
   * fixed part of every local file header, including in archives
   * created in a streaming fashion.
   */
  mtime: Date;
```

- [ ] **Step 4: Read the local-header dword in `Unzip.push`**

In `Unzip.push`, find this line (`src/index.ts:3781`):

```ts
            let lsc = b4(buf, i + 18), lsu = b4(buf, i + 22);
```

Add a line directly after it:

```ts
            // +10 is the DOS mod-time dword in a LOCAL header (the
            // central directory puts it at +12). Read it here, inside
            // the `l > i + 30 + fnl + es` guard, for two reasons:
            // outside the guard those 30 bytes are not yet known to be
            // buffered and b4 past the end silently yields zeros, and
            // the next line mutates i past the filename.
            const mtm = b4(buf, i + 10);
```

As in Tasks 1 and 2, `mtm` holds the RAW dword. `dosdt` is applied at
the point of use in the next step, so the name means the same thing in
all three call sites.

- [ ] **Step 5: Set it on the file object**

In the same block, extend the `file` object literal
(`src/index.ts:3788`). Before:

```ts
            const file = {
              name: fn,
              compression: cmp,
```

After:

```ts
            const file = {
              name: fn,
              compression: cmp,
              mtime: dosdt(mtm),
```

- [ ] **Step 6: Run the test to verify it passes**

```bash
npm run build && npm run test:browser
```

Expected: PASS. The new test name appears under both `plain` and `min`
with both assertions `ok`, and Tasks 1 and 2 still pass.

- [ ] **Step 7: Typecheck**

```bash
npx tsc --noEmit -p tsconfig.json > /dev/null; echo "EXIT:$?"
```

Expected: `EXIT:0`.

- [ ] **Step 8: Run the full suite**

```bash
npm test
```

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add src/index.ts test/browser/zip.ts
git commit -m "feat: expose entry mtime on the streaming Unzip class

Reads the DOS mod-time dword from each local file header at +10 and
sets it as UnzipFile.mtime."
```

---

## Task 4: Documentation and issue closeout

**Files:**
- Modify: `README.md` (:399-407 and :461-481)
- Modify: `TODO/README.md:12`

**Interfaces:**
- Consumes: the `mtime` members from Tasks 1 and 3.
- Produces: nothing.

There is no automated test for this task. Its verification is that the
code it documents already passes Tasks 1 to 3, plus a read-through.

- [ ] **Step 1: Update the `unzipSync` filter example**

In `README.md`, replace the filter block at :403-406. Before:

```js
  filter(file) {
    // Don't decompress the massive image or any files larger than 10 MiB
    return file.name != 'massiveImage.bmp' && file.originalSize <= 10_000_000;
  }
```

After:

```js
  filter(file) {
    // file.mtime is the entry's last modification time, decoded from the
    // ZIP header. DOS timestamps have two-second resolution and carry no
    // time zone, so it is a local-time Date on an even second.
    console.log(file.name, 'last modified', file.mtime);

    // Don't decompress the massive image or any files larger than 10 MiB
    return file.name != 'massiveImage.bmp' && file.originalSize <= 10_000_000;
  }
```

- [ ] **Step 2: Update the streaming `Unzip` example**

In `README.md`, replace the metadata lines at :471-474. Before:

```js
    // File sizes are sometimes not set if the ZIP file did not encode
    // them, so you may want to check that file.size != undefined
    console.log('Compressed size', file.size);
    console.log('Decompressed size', file.originalSize);
```

After:

```js
    // File sizes are sometimes not set if the ZIP file did not encode
    // them, so you may want to check that file.size != undefined
    console.log('Compressed size', file.size);
    console.log('Decompressed size', file.originalSize);

    // Unlike the sizes, mtime is always present: every local file
    // header carries a modification time
    console.log('Last modified', file.mtime);
```

- [ ] **Step 3: Tick the issue checkbox**

In `TODO/README.md`, line 12. Before:

```markdown
* [ ] [Expose entry mtime in unzip API](./mtime.md)
```

After:

```markdown
* [x] [Expose entry mtime in unzip API](./mtime.md)
```

- [ ] **Step 4: Verify the README edits landed where intended**

```bash
git diff --stat README.md TODO/README.md
grep -n "mtime" README.md
```

Expected: `README.md` and `TODO/README.md` both show as modified, and
`grep` reports the three new `mtime` mentions plus nothing unexpected.

- [ ] **Step 5: Run the full suite one final time**

```bash
npm test
```

Expected: PASS. Nothing in this task changes code, so a failure here
means an earlier task regressed.

- [ ] **Step 6: Commit**

```bash
git add README.md TODO/README.md
git commit -m "docs: document unzip entry mtime and close the issue"
```

---

## Out of scope

Do not do any of the following, all of which the spec explicitly ruled
out:

- Exporting `dosdt` as public API.
- Exposing the raw dword alongside the `Date`.
- Clamping or nulling an all-zero mod-time field. It decodes literally
  to 1979-11-30T00:00:00 local, and that is intended.
- Adding `mtime` to `Unzipped` (the `unzipSync` return value).
- Regenerating `docs/` typedoc output.
- Adding ZIP tests to `test/0-valid.ts`.
