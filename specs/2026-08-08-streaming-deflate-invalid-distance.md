# Streaming Deflate emits a distance into its own zero padding

Status: FIXED (US-013). Diagnosed under US-012, which shipped no source
or test changes; the "Proposed fix" section below was then applied
verbatim, with `test/10-stream-window.ts` as the regression guard.
Source: `TODO/streaming.md`

## Symptom

`Deflate` (and therefore `ZipDeflate`, `Gzip`, `Zlib`) can emit a
back-reference whose distance is larger than the number of bytes emitted
so far, which violates RFC 1951 section 3.2.5. Both fflate's own
`inflateSync` and node's `zlib.inflateRawSync` reject it.

Reproduced with the 234-byte EMF prefix from `TODO/streaming.md`, run
from a real `.mjs` file against `dist/node/index.js`:

```
--- level 3
  sync:   153B  fflate:OK  zlib:OK
  oneshot:153B  identical:true  zlib:OK
  stream: 153B  identical:false  fflate:FAIL: invalid distance
                                 zlib:FAIL: invalid distance too far back
  first divergence at byte 0: sync=0x4d stream=0x63
```

Levels 1, 3, 6 and 9 all reproduce. Level 1 also differs in length
(155 sync vs 156 streaming); levels 3, 6 and 9 match in length. So the
reported shape holds: `deflateSync` and callback `deflate()` are
byte-identical and valid, streaming is not.

## Responsible code

`dflt`, `src/index.ts:741`:

```ts
const maxd = Math.min(32767, i)
```

used as the loop guard at `src/index.ts:745`:

```ts
while (dif <= maxd && --ch && imod !== pimod) {
```

The state that makes that bound wrong is set up in `dfltInit`,
`src/index.ts:1415-1434`, and torn down in `Deflate.push`,
`src/index.ts:1495-1496`.

## Mechanism

`dflt` treats `i` as "number of bytes of real data that precede this
position", so `Math.min(32767, i)` is meant to say "never look further
back than the start of the stream". That is true for `deflateSync`,
where `dat` is the caller's buffer and index 0 is the first byte of the
stream.

It is false for the streaming classes. `dfltInit` allocates a
98304-byte scratch buffer and starts the cursor a full window in:

```ts
const s:DeflateState = { l:0, i:32768, w:32768, z:32768 }
const b = new U8(98304)   // 32768 lookback + 65536 chunk
```

Pushed data lands at `b[32768]` onward. `b[0..32767]` is reserved
lookback space that is **still zero-filled** until the buffer wraps.
`i` is therefore an index into that scratch buffer, and it overstates
the amount of emitted data by 32768. `maxd` lets the match finder walk
into the zero prefix, and a match found there is encoded as a perfectly
ordinary length/distance pair that no decompressor can resolve.

`deflateSync` gets this right for free: base index 0, so `i` really is
the byte count.

### The offending match

Instrumenting the `if (d)` branch of `dflt` to log every `(i, l, d)`,
sync and streaming agree on all 31 matches and disagree only on the
last one:

```
SYNC     ... M 218 3 72     M 229 4 184
STREAM   ... M 32986 3 72   M 32997 5 230
```

Subtract the 32768 base from the streaming indices and they line up
exactly: 218 -> 32986, 229 -> 32997. At relative offset 229 the sync
encoder picks length 4 at distance 184; the streaming encoder finds a
length-5 match at distance 230, which points at absolute `b[32767]` --
one byte before the data starts.

It is a genuine match *within the buffer*: `b[32767..32771]` is
`00 01 00 00 00` (the zero pad plus the input's first four bytes) and
the input's last five bytes are `00 01 00 00 00`. It is not a match
within the *stream*.

Instrumenting `inflt`'s dictionary branch confirms the decoder's view:

```
INFLBAD { bt: 229, dt: 230, add: 5, dl: 0 }
```

229 bytes produced, distance 230 requested.

## Why the divergence looks like it starts at byte 0

`TODO/streaming.md` reads the byte-0 difference as "a fundamental
divergence in the encoding decisions". It is not. The whole 234-byte
input is a single dynamic-Huffman block, so changing one symbol at the
very end changes the literal/length and distance frequencies, which
changes the Huffman code lengths, which changes the code-length tree
written in the block header. One bad symbol at relative offset 229
rewrites the block from byte 0. The 31 preceding matches are identical.

## Exposure window

The zero prefix only exists until the scratch buffer wraps.
`Deflate.push` (`src/index.ts:1489-1496`) does:

```ts
this.b.set(this.b.subarray(-32768))   // real data now fills b[0..32767]
this.b.set(chunk.subarray(split), 32768)
this.s.i = 32766
this.s.w = 32768
```

After that copy `b[0..32767]` holds genuinely emitted bytes and every
distance up to 32767 is legal. The bug is therefore live from the first
push until the cumulative pushed size exceeds `98304 - 32768 = 65536`
bytes, and dormant afterwards.

Empirically, over 300 random inputs of 1..5000 bytes at levels 1, 6 and
9, every failure came from the zero-heavy generator (about 60 percent
zero bytes) and none from the uniform-random one -- which is what you
would expect when the false match has to match a run of zeros. Larger
inputs (70000, 120000, 300000 bytes of random data) all passed. This
also explains the reporter's table: their all-zeros and repetitive-text
probes passed because a match against the zero prefix at distance
`d > i` is not *shorter* than the legal match the encoder would
otherwise pick, so it only wins on inputs where the tail happens to
align with the pad.

## Not the cause, but worth recording

`dopt` (`src/index.ts:1109`) derives the hash-table size from
`dat.length`. For `deflateSync` on this input that is 234 bytes, giving
`plvl = 12` (mask 4095). For the streaming path `dat` is the whole
98304-byte scratch buffer, giving `plvl = 18` (mask 262143, a 512 KB
`head` array per stream). Different hash tables can produce different
match choices, but here they did not: all 31 shared matches are
identical. This is a memory-use wart, not the correctness bug.

## Relationship to US-010

None. US-010 was a bit-cursor unpacking mistake in `Deflate.flush`'s
sync-flush marker, confined to `wfblk`'s `pos` argument. This one is
LZ77 window bookkeeping in `dflt`'s match finder. The two share no code
and the US-011 fix does not move this bug -- the failure reproduces
unchanged on the current tree.

## Proposed fix (for US-013)

`DeflateState` needs to carry the index of the first valid byte in
`dat` so `dflt` can bound the search by real bytes rather than by
buffer offset. `dfltInit` is the single construction point for all
three streaming classes (`Deflate` at `src/index.ts:1452`, `Gzip` at
`1879`, `Zlib` at `2245`), so one field covers them all.

```ts
// DeflateState
b?:number;              // index of first valid byte in dat

// dflt, src/index.ts:741
const maxd = Math.min(32767, i - (st.b || 0))

// dfltInit, src/index.ts:1424
const s:DeflateState = { l:0, i:32768, w:32768, z:32768, b:32768 }
...
s.i = s.b = 32768 - dict.length   // dictionary bytes ARE legal targets

// Deflate.push, after the wrap copy at src/index.ts:1495
this.s.b = 0
```

`|| 0` keeps `deflateSync` and the one-shot paths at base 0 with no
other change. The dictionary case sets the base to the start of the
dictionary rather than to 32768, because a decompressor primed with the
same dictionary can resolve those distances -- this matches what `dopt`
already does for the non-streaming dictionary path, where the
dictionary is prepended to `dat` and so sits at base 0.

Note the name collision hazard: `InflateState` also has a `b` member
(`src/index.ts:276`, the output byte count). They are different types
and never mix, but a future reader will trip over it. `Deflate.flush`
does not need a change -- it only calls `this.p`, and the base is
unaffected by a flush.

### Verification of the proposed fix

The change above was applied temporarily and then reverted, since
US-012 ships no source changes. With it applied:

- EMF prefix, levels 1-9 crossed with chunk shapes
  {one push, 1-byte, 7-byte, 64-byte, 200-byte}: 45/45 pass under both
  `zlib.inflateRawSync` and fflate's `inflateSync`. Without it, 45/45
  fail.
- 300 generated inputs (uniform random, 60 percent zeros, periodic) of
  1..5000 bytes at levels 1, 6, 9: 0 failures. Without it, 24 failures,
  all from the zero-heavy generator.
- 70000/120000/300000-byte random streams at one push and 4096-byte
  chunks: 0 failures both with and without.
- `deflateSync` output unchanged; `npx tsx test/0-valid.ts` 13/13.
