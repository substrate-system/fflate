import { inflateRawSync } from 'node:zlib'
import { test } from '@substrate-system/tapzero'
import type { Test } from '@substrate-system/tapzero'
import {
    Deflate,
    type DeflateOptions,
    Zip,
    ZipDeflate,
    deflateSync,
    inflateSync,
    unzipSync
} from 'fflate'

// Regression coverage for TODO/streaming.md, diagnosed in
// specs/2026-08-08-streaming-deflate-invalid-distance.md.
//
// The streaming classes deflate out of a 98304-byte scratch buffer whose
// first 32768 bytes are reserved lookback space, zero-filled until the
// buffer wraps. `dflt` bounded the match search with
// `Math.min(32767, i)`, treating the scratch-buffer index as a count of
// emitted bytes, so the match finder was allowed to walk into that zero
// pad. A match found there encodes as an ordinary length/distance pair
// whose distance exceeds the bytes emitted so far -- an RFC 1951
// section 3.2.5 violation that no inflater can resolve.
//
// The input below is the 234-byte EMF prefix from the bug report. It is
// the specific case where the false match against the pad beats the
// legal one: its last five bytes are `00 01 00 00 00`, which matches the
// final pad byte followed by the input's first four bytes.

const EMF_PREFIX_B64 =
    'AQAAAGwAAAAAAAAAAAAAAM0AAAAcAAAAAAAAAAAAAABqHAAA/gMAACBFTUYAAAEAAAsAADIAAAAH' +
    'AAAAAAAAAAAAAAAAAAAAAAUAAAAEAADEAQAAaQEAAAAAAAAAAAAAAAAAAOPjBgAcgwUARgAAAGQC' +
    'AABWAgAAR0RJQwEAAIAAAwAAfQ+OiAAAAAA+AgAAAQAJAAADHwEAAAYAKwAAAAAABAAAAAMBCAAF' +
    'AAAACwIAAAAABQAAAAwCHQDOAAMAAAAeAAcAAAD8AgAAaWlpAAAABAAAAC0BAAAJAAAAHQYhAPAA' +
    'HQABAAAA'

const EMF = Uint8Array.from(Buffer.from(EMF_PREFIX_B64, 'base64'))

// The bug is level- and chunk-shape independent per the report, so the
// matrix below crosses every non-stored level it names with the two
// extreme chunk shapes: everything in one push, and a push per byte.
type Level = DeflateOptions['level']
const LEVELS:Level[] = [1, 3, 6, 9]

const join = (chunks:Uint8Array[]) => {
    let n = 0
    for (const c of chunks) n += c.length
    const out = new Uint8Array(n)
    let i = 0
    for (const c of chunks) {
        out.set(c, i)
        i += c.length
    }
    return out
}

// Split `input` into chunks of `size` bytes; size 0 means one chunk.
const chunkUp = (input:Uint8Array, size:number) => {
    if (!size) return [input]
    const out:Uint8Array[] = []
    for (let i = 0; i < input.length; i += size) {
        out.push(input.subarray(i, Math.min(i + size, input.length)))
    }
    return out
}

const streamDeflate = (input:Uint8Array, level:Level, chunk:number) => {
    const chunks:Uint8Array[] = []
    const d = new Deflate({ level }, (c) => {
        if (c.length) chunks.push(new Uint8Array(c))
    })
    const parts = chunkUp(input, chunk)
    for (let i = 0; i < parts.length; ++i) {
        d.push(parts[i], i === parts.length - 1)
    }
    return join(chunks)
}

// Both inflaters have to agree, because fflate's own decoder rejecting
// the stream would not prove anything about interoperability, and node's
// rejecting it would not prove fflate can still read its own output.
const roundTrips = (
    t:Test,
    out:Uint8Array,
    expect:Uint8Array,
    why:string
) => {
    t.deepEqual(
        Array.from(inflateSync(out)),
        Array.from(expect),
        `${why}: fflate inflateSync recovers the input`
    )
    t.deepEqual(
        Array.from(new Uint8Array(inflateRawSync(out))),
        Array.from(expect),
        `${why}: zlib inflateRawSync recovers the input`
    )
}

test('streaming Deflate never emits an out-of-window distance', t => {
    for (const level of LEVELS) {
        for (const chunk of [0, 1]) {
            const shape = chunk ? '1-byte pushes' : 'one push'
            roundTrips(
                t,
                streamDeflate(EMF, level, chunk),
                EMF,
                `EMF prefix, level ${level}, ${shape}`
            )
        }
    }
})

// The bug was originally found through ZipDeflate, which wraps Deflate
// and so inherits its state. Driving it through a real Zip keeps the
// coverage on the path the reporter actually used.
test('ZipDeflate output for the same input is a readable archive', t => {
    for (const level of LEVELS) {
        const chunks:Uint8Array[] = []
        const z = new Zip((e, d) => {
            if (e) throw e
            if (d.length) chunks.push(new Uint8Array(d))
        })
        const f = new ZipDeflate('image2.emf', { level })
        z.add(f)
        f.push(EMF, true)
        z.end()

        const files = unzipSync(join(chunks))
        t.deepEqual(
            Array.from(files['image2.emf']),
            Array.from(EMF),
            `ZipDeflate level ${level}: unzipSync recovers the input`
        )
    }
})

// The fix bounds the search by the first valid byte of the scratch
// buffer, which is index 0 for every non-streaming path. Nothing about
// deflateSync's output should move. test/1-size.ts holds the byte-size
// expectations for the real corpus; this asserts the invariant directly
// on the input that exercised the bug.
test('deflateSync output is unchanged by the window fix', t => {
    // Measured on the tree before the fix, in level order 1, 3, 6, 9.
    const expected = [155, 153, 153, 153]
    LEVELS.forEach((level, i) => {
        const out = deflateSync(EMF, { level })
        t.equal(
            out.length,
            expected[i],
            `deflateSync level ${level} still emits ${expected[i]} bytes`
        )
        roundTrips(t, out, EMF, `deflateSync level ${level}`)
    })
})

// A window base that is too SMALL would silently cost compression
// instead of corrupting output, so assert the streaming encoder still
// reaches back across a chunk boundary. 40000 bytes forces the scratch
// buffer to wrap, which is where the base resets to 0.
test('streaming Deflate still finds matches across the window', t => {
    const unit = new Uint8Array(1000)
    for (let i = 0; i < unit.length; ++i) unit[i] = (i * 7 + 3) & 255
    const input = new Uint8Array(80000)
    for (let i = 0; i < 80; ++i) input.set(unit, i * 1000)

    const out = streamDeflate(input, 6, 4096)
    roundTrips(t, out, input, 'repeating 1 KB unit over 80 KB')
    t.ok(
        out.length < 4000,
        `long-range matches still found (${out.length} bytes for 80 KB)`
    )
})
