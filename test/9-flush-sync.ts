import { inflateRawSync } from 'node:zlib'
import { test } from '@substrate-system/tapzero'
import type { Test } from '@substrate-system/tapzero'
import { Deflate, inflateSync, strToU8, strFromU8 } from 'fflate'

// Regression coverage for TODO/flush.md, diagnosed in
// specs/2026-08-08-flush-sync-corruption.md.
//
// `Deflate.flush(true)` appends an empty non-final stored block as a
// sync marker. It passed the PACKED `st.r` -- bit offset in bits 0-2,
// leftover partial byte in bits 3+ -- to `wfblk`, which wants a plain
// bit position, and it dropped the `+ 1` for the BFINAL bit that every
// other `wfblk` call site accounts for. The marker came out with
// LEN/NLEN at the wrong offset, or truncated entirely, and inflaters
// rejected it with "invalid stored block lengths".
//
// The failure is bit-alignment dependent, not size dependent: it fires
// only when the pending bit offset is 6, or when the leftover partial
// byte is non-zero. The sizes below are the ones the diagnosis recorded
// as landing on each of those cases.

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

// Push `input` in one chunk, flush, then terminate the stream.
const compressWithFlush = (input:Uint8Array, sync:boolean) => {
    const chunks:Uint8Array[] = []
    const d = new Deflate((c) => chunks.push(c))
    d.push(input)
    d.flush(sync)
    d.push(new Uint8Array(0), true)
    return join(chunks)
}

// Push each chunk followed by its own sync flush.
const compressMultiFlush = (parts:Uint8Array[]) => {
    const chunks:Uint8Array[] = []
    const d = new Deflate((c) => chunks.push(c))
    for (const p of parts) {
        d.push(p)
        d.flush(true)
    }
    d.push(new Uint8Array(0), true)
    return join(chunks)
}

const roundTrips = (t:Test, out:Uint8Array, expect:Uint8Array, why:string) => {
    t.deepEqual(
        Array.from(inflateSync(out)),
        Array.from(expect),
        `${why}: fflate inflateSync recovers the input`
    )
    t.deepEqual(
        Array.from(new Uint8Array(inflateRawSync(Buffer.from(out)))),
        Array.from(expect),
        `${why}: zlib inflateRawSync recovers the input`
    )
}

// n = 62472 is the smallest failing size the diagnosis found (bit
// offset exactly 6, zero partial byte). 70000 and 200000 leave a
// non-zero partial byte. 1000 and 62471 pass even before the fix, and
// are here to prove the fix does not regress the aligned cases.
const SIZES = [0, 1, 1000, 62471, 62472, 70000, 100000, 200000]

for (const n of SIZES) {
    test(`flush(true) on ${n} repeated bytes stays inflatable`, (t:Test) => {
        const input = strToU8('a'.repeat(n))
        roundTrips(t, compressWithFlush(input, true), input, `n=${n} sync`)
    })
}

test('flush(false) still round-trips', (t:Test) => {
    for (const n of SIZES) {
        const input = strToU8('a'.repeat(n))
        roundTrips(t, compressWithFlush(input, false), input, `n=${n} async`)
    }
})

test('repeated sync flushes stay inflatable', (t:Test) => {
    // Each flush lands on a different bit alignment, so a sequence
    // exercises several of them in one stream.
    const parts = [
        strToU8('a'.repeat(70000)),
        strToU8('b'.repeat(9)),
        strToU8('hello world'),
        strToU8('c'.repeat(40000))
    ]
    const expect = join(parts)
    roundTrips(t, compressMultiFlush(parts), expect, 'multi-flush')
    t.equal(
        strFromU8(inflateSync(compressMultiFlush(parts))),
        strFromU8(expect),
        'multi-flush: text survives intact'
    )
})

test('sync flush markers scan clean across many alignments', (t:Test) => {
    // The bug is a bit-alignment coincidence, so a sweep is the only
    // honest check. Mixed bytes shift the alignment differently than
    // the runs above do.
    let bad = 0
    for (let n = 0; n < 600; ++n) {
        const input = seqBytes(n)
        const out = compressWithFlush(input, true)
        try {
            const got = new Uint8Array(inflateRawSync(Buffer.from(out)))
            if (got.length !== n) bad++
            else {
                for (let i = 0; i < n; ++i) {
                    if (got[i] !== input[i]) { bad++; break }
                }
            }
        } catch (_err) {
            bad++
        }
    }
    t.equal(bad, 0, 'every size in 0..599 round-trips through zlib')
})

// Pseudo-random but deterministic bytes, so alignments vary by size.
function seqBytes (n:number):Uint8Array {
    const out = new Uint8Array(n)
    let s = 0x2545f491
    for (let i = 0; i < n; ++i) {
        s ^= s << 13; s >>>= 0
        s ^= s >>> 17
        s ^= s << 5; s >>>= 0
        out[i] = s & 255
    }
    return out
}
