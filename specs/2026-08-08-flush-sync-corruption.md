# flush(true) produces an invalid empty stored block

Status: FIXED in US-011. Diagnosed in US-010; this document is kept as the
record of the mechanism. Regression coverage lives in
`test/9-flush-sync.ts`.
Source: `TODO/flush.md`

## Symptom

```js
const chunks = []
const d = new Deflate((c) => chunks.push(c))
d.push(new TextEncoder().encode('a'.repeat(70000)))
d.flush(true)
d.push(new Uint8Array(0), true)
zlib.inflateRawSync(Buffer.concat(chunks))  // invalid stored block lengths
```

`flush(false)` on the same input round-trips fine, so the compressed
blocks themselves are correct. Only the extra sync-flush marker that
`flush(true)` appends is malformed.

## Responsible code

`Deflate.prototype.flush`, `src/index.ts:1522-1538`, specifically the
`sync` branch at `src/index.ts:1530-1537`:

```ts
if (sync) {
    const c = new U8(6)
    c[0] = this.s.r >> 3
    // write empty, non-final type-0 block
    const ep = wfblk(c, this.s.r, et)
    this.s.r = 0
    this.ondata(c.subarray(0, ep >> 3), false)
}
```

The supporting pieces are `wfblk` (`src/index.ts:581-591`) and the place
`st.r` is produced, `dflt` (`src/index.ts:797-806`).

## What `st.r` actually holds

At the end of a non-final `dflt` call the bit cursor `pos` usually sits
mid-byte. `dflt` cannot emit that partial byte, so it stashes it:

```ts
st.r = (pos & 7) | w[(pos / 8) | 0] << 3
pos -= 7
```

`st.r` is therefore a PACKED pair, not a bit position:

- bits 0-2 -- the bit offset within the partial byte, `pos & 7`
- bits 3+  -- the value of the partial byte itself

`flush` unpacks the byte half correctly (`c[0] = this.s.r >> 3`) and then
passes the whole packed value straight into `wfblk` as if it were a bit
position. That is the defect.

## What `wfblk` expects

```ts
const wfblk = (out, pos, dat) => {
    const s = dat.length
    const o = shft(pos + 2)
    out[o]     = s & 255
    out[o + 1] = s >> 8
    out[o + 2] = out[o] ^ 255
    out[o + 3] = out[o + 1] ^ 255
    ...
    return (o + 4 + s) * 8
}
```

Its `pos` is the bit position immediately AFTER the BFINAL bit. Every
other call site honours that: `dflt` writes BFINAL itself and then calls
`wfblk(w, pos + 1, ...)` (`src/index.ts:815`). `pos + 2` accounts for the
two BTYPE bits, and `shft` rounds up to the byte boundary where LEN
begins.

`flush` never adds the `+ 1`, and it passes the packed `st.r` rather than
`st.r & 7`.

## The two resulting failure modes

Instrumented run, printing the packed `st.r` at the flush point and the
bytes the marker actually emits:

```
      n   st.r  bitpos  byte   marker bytes      inflateRawSync
      1      2       2     0   00 00 00 ff ff    OK
   1000      3       3     0   00 00 00 ff ff    OK
  40000      5       5     0   00 00 00 ff ff    OK
  62471      5       5     0   00 00 00 ff ff    OK
  62472      6       6     0   00 00 00 ff ff    invalid stored block lengths
  63472      3       3     0   00 00 00 ff ff    OK
  70000     19       3     2   02 00 00 00 00 ff invalid stored block lengths
 100000      6       6     0   00 00 00 ff ff    invalid stored block lengths
 200000     19       3     2   02 00 00 00 00 ff invalid stored block lengths
```

### Mode 1 -- the leftover bits are non-zero (`st.r >> 3 != 0`)

`n = 70000`: `st.r = 19`, so the real bit offset is 3 and the partial
byte is `0x02`. `wfblk` receives `pos = 19` and computes
`o = shft(21) = 3`, so it writes LEN at `c[3]`/`c[4]` and NLEN at
`c[5]`/`c[6]` -- but `c` is only 6 bytes, so `c[6]` is dropped silently
(a typed-array out-of-range write is a no-op, no throw). The marker
becomes `02 00 00 00 00 ff`: LEN reads as `0x0000` and NLEN as `0x00ff`.
`NLEN != ~LEN`, so the inflater rejects it.

Whenever the partial byte is non-zero the marker is wrong. `o` scales
with the byte value, so for larger partial bytes `o` exceeds 6 entirely
and every LEN/NLEN write is discarded, leaving five or six zero bytes.

### Mode 2 -- bit offset exactly 6, even with a zero partial byte

`n = 62472` and `n = 100000`: `st.r = 6`, so the byte value is 0 and the
packed value coincidentally equals the bit offset. `wfblk` still drops
the `+ 1`, computing `o = shft(6 + 2) = shft(8) = 1` instead of
`shft(6 + 3) = shft(9) = 2`.

Six bits are already consumed in `c[0]`, so the 3-bit block header
occupies bits 6 and 7 of `c[0]` plus bit 0 of `c[1]`. LEN must therefore
start at `c[2]`, not `c[1]`. The marker is emitted one byte short and the
inflater reads LEN/NLEN from the wrong offset.

Offset 6 is the only value where the missing `+ 1` matters on its own:
`shft(p + 2)` and `shft(p + 3)` agree for `p` in 0..5 and for `p = 7`,
and differ only at `p = 6`.

## Why the bug looks size-dependent

It is not a size threshold, it is a data-dependent bit-alignment
coincidence. The output is correct exactly when the bit offset is in
`{0, 1, 2, 3, 4, 5, 7}` AND the leftover partial byte is zero. Small
inputs happen to land there.

A binary search over `'a'.repeat(n)` gives a first failure at
**n = 62472** (62471 passes), but that is not a threshold: 63472 passes
again while 70000, 100000 and 200000 all fail. Recording it because
US-010 asks for it, and flagging that US-011's "size below the
threshold" case should be understood as "a size that happens to pass",
not as a safe region.

## Correct bytes per RFC 1951 section 3.2.4

A non-final empty stored block is:

1. BFINAL, 1 bit, value 0.
2. BTYPE, 2 bits, value 00.
3. Skip any remaining bits to the next byte boundary.
4. LEN, 2 bytes little-endian, value `0x0000`.
5. NLEN, 2 bytes little-endian, the one's complement of LEN, `0xFFFF`.
6. Zero data bytes.

So on an already byte-aligned stream the marker is the canonical
`00 00 00 FF FF`. With `k` bits already pending in the partial byte, the
partial byte is emitted first with the three header bits ORed in at bit
`k` (all zero, so it is unchanged), then padding zeros to the boundary,
then `00 00 FF FF`. For `k` in 1..5 that is 5 bytes total; for `k` of 6
or 7 the header spills into the next byte and it is 6 bytes.

## Fix

Pass the unpacked bit offset and restore the BFINAL bit that every other
`wfblk` call site accounts for:

```ts
const ep = wfblk(c, (this.s.r & 7) + 1, et)
```

BFINAL is 0 here so it does not need to be written, only counted.
`new U8(6)` is already large enough: the worst case is offset 6 or 7,
giving `o = 2` plus 4 bytes of LEN/NLEN.

Verified by applying that single-line change and sweeping
`'a'.repeat(n)` for n in 0..2000 plus every 617th size up to 300000
(about 2500 sizes), asserting both `zlib.inflateRawSync` and fflate's
own `inflateSync` recover the exact input. Zero failures. The change was
reverted afterwards -- US-010 ships no source changes.

## Relationship to US-012

Probably none. This is purely a bit-cursor unpacking mistake in the
sync-flush marker; it does not touch the LZ77 window, the hash chains,
or `st.i` / `st.w`. US-012's out-of-window distance is a window
bookkeeping problem and should be diagnosed independently.
