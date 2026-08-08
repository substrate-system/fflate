# Expose entry mtime in the unzip API

Date: 2026-08-08
Status: approved, not yet implemented

Issue: [`TODO/mtime.md`](../TODO/mtime.md), from
[101arrowz/fflate#290](https://github.com/101arrowz/fflate/issues/290).

## Goal

Every ZIP header already carries a last-modification timestamp, and
`fflate` already *writes* one on the zip side via `ZipAttributes.mtime`.
The unzip side discards it. Consumers who want it are reduced to patching
the library in their bundler and re-deriving the bit layout by hand.

Expose it as a decoded `Date` at the three points where unzip surfaces
per-entry metadata:

- the `filter` callback in `unzipSync` (`UnzipFileInfo`)
- the `filter` callback in async `unzip` (`UnzipFileInfo`)
- the file objects emitted by the streaming `Unzip` class (`UnzipFile`)

## The format

DOS date-time is a single 32-bit little-endian word, packed as:

```
 31        25 24    21 20   16 15    11 10       5 4        0
+------------+--------+-------+--------+----------+----------+
| year-1980  | month  |  day  |  hour  |  minute  |  sec / 2 |
|   7 bits   | 4 bits | 5 bits| 5 bits |  6 bits  |  5 bits  |
+------------+--------+-------+--------+----------+----------+
```

Two consequences follow from the layout and both shape the design:

1. **Two-second resolution, no milliseconds.** The seconds field stores
   `seconds / 2`. A round-trip through a ZIP header cannot preserve odd
   seconds or sub-second precision.
2. **No time zone.** The value is bare wall-clock. `wzh`
   (`src/index.ts:2917`) writes it with local-time getters
   (`getFullYear`, `getMonth`, `getDate`, `getHours`, `getMinutes`,
   `getSeconds`), so it must be read back with the local-time `Date`
   constructor for a `zipSync` -> `unzipSync` round-trip to be exact.
   Decoding as UTC would shift every timestamp by the runner's offset.

## Design

### 1. The decode helper

A single internal, unexported helper, placed next to `zh`/`z64hs` around
`src/index.ts:2862`:

```ts
// decode a DOS date-time dword into a Date
const dosdt = (t: number) => new Date(
  ((t >>> 25) & 127) + 1980,
  ((t >>> 21) & 15) - 1,
  (t >>> 16) & 31,
  (t >>> 11) & 31,
  (t >>> 5) & 63,
  (t & 31) * 2
);
```

Two details are load-bearing.

**`>>>`, not `>>`.** The year field is 7 bits wide at bit 25, so any
year from 2044 on sets bit 31 -- 2099 packs as `0xEF9FBF7D`. Under `>>`
the value is coerced to a signed int32 and the shift sign-extends.

Measured: with the `& 127` in place, `>>` and `>>>` both yield 2099, so
the mask alone is sufficient *as written*. The unsigned shift is not
therefore redundant. `((t >>> 25) & 127) + 1980` invites simplification
to `(t >>> 25) + 1980`, which is correct only with the unsigned shift;
the signed version of that same simplification returns 1971. Keeping
the shift unsigned means the mask is a redundancy rather than the sole
thing holding the line. `b4` already returns unsigned via `>>> 0`, so
`dosdt` stays unsigned end to end.

No test distinguishes `>>` from `>>>` here, and none can while the mask
is present. That is a documented property of the decoder, not a gap in
the suite.

**Mask before subtracting.** The snippet in the original issue reads
`(time >> 21) & 0x0f - 1`, which parses as `& (0x0f - 1)` -- `& 14` --
because `-` binds tighter than `&` in JavaScript. That silently clears the
low bit of the month. The correct form parenthesises the mask, as above.

It takes the raw dword rather than a `(d, b)` pair because the three call
sites read from different headers at different offsets.

### 2. Reading the raw dword

The field sits at **+12** in a central directory header and at **+10** in
a local file header. Both offsets follow from `wzh`, which writes a
central header as `sig(4) versionMadeBy(2) versionNeeded(2) flags(2)
compression(2) time(4)` and a local header as `sig(4) version(2) flags(2)
compression(2) time(4)`.

**`zh()` (`src/index.ts:2865`)** gains a 7th element in its returned
tuple: `b4(d, b + 12)` -- the **raw** dword, not a `Date`.

Keeping it raw is deliberate. Both `unzip` and `unzipSync` guard the
filter with `!fltr || fltr({...})`, so short-circuit evaluation means the
object literal is never constructed when no filter is supplied. Returning
a `Date` from `zh` instead would allocate one per central directory entry
on every unzip, filter or not.

**`Unzip.push` (`src/index.ts:3788`)** reads the local header directly
rather than going through `zh`. The read goes *inside* the
`if (l > i + 30 + fnl + es)` guard, alongside `lsc`/`lsu`, for two
reasons: outside the guard the 30 header bytes are not yet known to be
buffered (`b4` past the end silently yields zeros), and `i` is mutated by
the filename slice on the following line, so the read must precede it.

### 3. Zip64 and data descriptors

Neither needs special handling.

Zip64 relocates only the compressed size, uncompressed size, and local
header offset into the extra field -- that is exactly the set `z64hs`
handles, and the time field is not among them. An entry with the data
descriptor flag (bit 3) defers only its CRC and sizes to a trailing
descriptor record; the time is still written in the local header at the
usual offset. In both cases `+12`/`+10` remain correct.

### 4. Public API

`UnzipFileInfo` and `UnzipFile` each gain the same member:

```ts
  /**
   * The last modification time of the file. DOS timestamps carry no
   * time zone and have two-second resolution, so this is a local-time
   * Date on an even second, with no milliseconds.
   */
  mtime: Date;
```

**Non-optional on both.** On `UnzipFile`, `size` and `originalSize` are
optional because archives produced in a streaming fashion omit them from
the local header -- but the DOS time field is a fixed part of every local
header and is always readable. Making `mtime` optional would force a
narrowing check on every consumer for a case that cannot occur.

**An all-zero field decodes literally.** Some tools write zeros to mean
"no timestamp". The decode does not special-case it: month `0 - 1` and
day `0` both roll backward, yielding `1979-11-30T00:00:00` local. This
matches `@zip.js`, which also does not special-case it, and it keeps the
type non-optional and the hot path branch-free. Consumers that care can
compare against that sentinel.

**`Unzipped` is unchanged.** The value returned by `unzipSync` is a plain
`name -> Uint8Array` map. There is nowhere to hang metadata without
changing its shape, and the issue does not ask for it.

## Testing

All tests go in `test/browser/zip.ts`. That file is the home of every ZIP
test in the project (`test/0-valid.ts` has no ZIP coverage), and the
browser suite runs each suite twice, once against the plain bundle and
once against the minified one.

Every test date uses an even second and zero milliseconds, since the
format holds neither.

1. **Per-entry round-trip through `unzipSync`.** Build an archive with
   `zipSync` giving two entries *different* `mtime` values via the
   `[data, opts]` tuple form, then assert each comes back distinctly
   through the filter. Distinct values are what prove the read is
   per-entry rather than one header being read for all of them.
2. **Same assertion through async `unzip`.** Exercises the same `zh`
   path but the async control flow.
3. **Same assertion through the streaming `Unzip` class.** This reads a
   different header at a different offset and is therefore genuinely
   separate coverage, not a duplicate.
4. **Boundary values.** `1980-01-01T00:00:00` local -- the DOS epoch, all
   packed fields at their minimum -- and `2099-12-31T23:59:58` local.

   2099 is the real top of the range -- `wzh` throws error 10 for any
   year past it -- and it packs as `0xEF9FBF7D`, with bit 31 set. The
   case pins the year field's shift distance and mask width at the one
   value where an off-by-one in either is visible. It does **not**
   distinguish `>>` from `>>>` in `dosdt`; see the note under the
   decoder above for why nothing can.

## Out of scope

- Exporting `dosdt` as public API.
- Exposing the raw dword alongside the `Date`.
- Clamping or nulling an all-zero timestamp.
- Adding `mtime` to `Unzipped`.
- `docs/` is generated typedoc, and `npm run build-docs` is on the
  project's never-run list. It regenerates on the next docs build.

## Files touched

| File | Change |
| --- | --- |
| `src/index.ts` | `dosdt` helper; `zh` returns the raw dword; `mtime` on `UnzipFileInfo`, `UnzipFile`, and the three construction sites |
| `test/browser/zip.ts` | The four tests above |
| `README.md` | `mtime` in the existing filter example (~:403) and streaming example (~:474) |
| `TODO/README.md` | Tick the checkbox |
