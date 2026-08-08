**fflate version:** 0.8.3 (also reproduces on 0.8.2)
**Node version:** v22.22.0
**OS:** Linux

## The problem

`Deflate` (the low-level streaming class) and `ZipDeflate` (which wraps it) produce a **deflate stream that off-the-shelf inflaters reject as invalid** for some inputs. `deflateSync` and the callback-style `deflate()` on the same input are both correct and byte-identical.

The streaming-produced bytes diverge from the sync output starting **at offset 0** — the deflate block-type bits themselves differ, so it's a fundamental divergence in the streaming class's encoding decisions, not a chunk-boundary or window edge case. The result emits at least one back-reference whose distance is illegal per RFC 1951.

Discovered while building an XLSB → XLSX converter that streams ZIP output via `ZipDeflate`. One workbook produced a ZIP whose `xl/media/image2.emf` entry failed to extract. The trigger was a specific 234-byte prefix of that EMF, bisected from the original 2816 bytes — trimming any further (from either end) makes the bug stop reproducing.

## Reproducer (self-contained, no extra deps)

```js
import { Deflate, deflate, deflateSync, inflateSync } from 'fflate';
import zlib from 'node:zlib';

// 234 bytes — prefix of a real Windows Enhanced Metafile (.emf)
// embedded in an XLSB workbook. Bisected from the full 2816-byte
// EMF; trimming either end makes the bug stop triggering. Contains
// only the EMF format magic + standard GDI record opcodes.
const TEST_INPUT = Uint8Array.from(Buffer.from(
  'AQAAAGwAAAAAAAAAAAAAAM0AAAAcAAAAAAAAAAAAAABqHAAA/gMAACBFTUYAAAEAAAsAADIAAAAH' +
  'AAAAAAAAAAAAAAAAAAAAAAUAAAAEAADEAQAAaQEAAAAAAAAAAAAAAAAAAOPjBgAcgwUARgAAAGQC' +
  'AABWAgAAR0RJQwEAAIAAAwAAfQ+OiAAAAAA+AgAAAQAJAAADHwEAAAYAKwAAAAAABAAAAAMBCAAF' +
  'AAAACwIAAAAABQAAAAwCHQDOAAMAAAAeAAcAAAD8AgAAaWlpAAAABAAAAC0BAAAJAAAAHQYhAPAA' +
  'HQABAAAA', 'base64'));

function streamingDeflate(buf, level) {
  const chunks = [];
  const d = new Deflate({ level }, (c) => { if (c?.length) chunks.push(new Uint8Array(c)); });
  d.push(buf, true);
  let total = 0; for (const c of chunks) total += c.length;
  const out = new Uint8Array(total); let off = 0;
  for (const c of chunks) { out.set(c, off); off += c.length; }
  return out;
}

const LEVEL = 3;
const sync = deflateSync(TEST_INPUT, { level: LEVEL });
const oneshot = await new Promise((res, rej) =>
  deflate(TEST_INPUT, { level: LEVEL }, (e, o) => e ? rej(e) : res(o)));
const stream = streamingDeflate(TEST_INPUT, LEVEL);

console.log(`deflateSync: ${sync.length}B  -- fflate.inflateSync: ${tryInflate(inflateSync, sync)}  zlib.inflateRawSync: ${tryInflate(zlib.inflateRawSync, sync)}`);
console.log(`deflate():   ${oneshot.length}B  -- fflate.inflateSync: ${tryInflate(inflateSync, oneshot)}  zlib.inflateRawSync: ${tryInflate(zlib.inflateRawSync, oneshot)}  identical-to-sync: ${eq(oneshot, sync)}`);
console.log(`Deflate():   ${stream.length}B  -- fflate.inflateSync: ${tryInflate(inflateSync, stream)}  zlib.inflateRawSync: ${tryInflate(zlib.inflateRawSync, stream)}  identical-to-sync: ${eq(stream, sync)}`);

function tryInflate(fn, b) { try { fn(b); return 'OK'; } catch (e) { return `FAIL: ${e.message}`; } }
function eq(a, b) { return a.length === b.length && Buffer.compare(Buffer.from(a), Buffer.from(b)) === 0; }
```

## Output (fflate 0.8.3, Node 22.22.0)

```
deflateSync: 153B  -- fflate.inflateSync: OK  zlib.inflateRawSync: OK
deflate():   153B  -- fflate.inflateSync: OK  zlib.inflateRawSync: OK  identical-to-sync: true
Deflate():   153B  -- fflate.inflateSync: FAIL: invalid distance  zlib.inflateRawSync: FAIL: invalid distance too far back  identical-to-sync: false
```

The streaming-emitted distance is greater than the actual amount of data fed so far (RFC 1951 §3.2.5 violation), so two independent inflaters reject it. `Deflate({level: 0}, ...)` (no compression) is fine because the stored encoding skips the back-reference path entirely.

## Audit — ruling out "wrong usage"

I went through this carefully before filing. The bug reproduces independently of:

| Variable | Variants tried | Same invalid output? |
|---|---|---|
| Input construction | `Uint8Array.from(Buffer.from(b64))`; raw `Buffer.from(b64)`; fresh `ArrayBuffer`-backed copy; `readFileSync` round-trip | yes — all 4 produce the same invalid 153 B |
| Chunk shape | single `push(buf, true)`; 1-byte chunks; 64-byte chunks; 200-byte chunk straddling the input | yes — all produce the same invalid 153 B |
| Level | 1, 3, 6, 9 (all non-zero) | yes — all reproduce; level 0 is fine (stored) |
| Inflater | fflate's `inflateSync`; Node's native `zlib.inflateRawSync` | both reject |

And the two non-streaming API paths (`deflateSync`, callback-style `deflate()`) produce byte-identical correct output for the same input. So the streaming `Deflate` class is the only piece that produces wrong bytes.

## Characterization — what does and doesn't trigger

Across the same 234-byte input shape, plus a handful of probes:

| Input (~2.8 KB sample sizes) | Streaming bytes match `deflateSync`? | Streaming inflates? |
|---|---|---|
| Random binary (`crypto.randomBytes`) | no — streaming output is slightly larger (~1 % overhead) | **yes** (valid but less efficient) |
| All-zeros buffer | yes | yes |
| Highly repetitive text | yes | yes |
| **Real-world binary like this EMF** | yes (same size, different bytes) | **no** |

So the bug isn't size- or level-dependent, and isn't a chunking edge case. It only fires for certain byte distributions.

## Suspected area

Without spelunking the deflate source: the symptom looks like the LZ77 match-finder accepting a candidate whose distance is greater than the spec-allowed maximum *or* greater than the actual amount of data fed so far. `deflateSync` handles it correctly, so whatever bookkeeping `deflateSync` does to bound the search window is missing or off-by-something in the streaming `Deflate` path. The first-byte divergence suggests the block-type choice or initial dictionary setup is different too.

Happy to dig further if it helps.

## Workaround used downstream

Routing all "I have the full bytes already" entries through `deflateSync` (and only using streaming `Deflate` for genuinely streamed payloads, where the bug hasn't surfaced for our inputs) avoids it. That's what we shipped in our converter while waiting on a fix.

Thanks for the library — `fflate` is great for the 99 % of cases that don't hit this.
