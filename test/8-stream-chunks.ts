import { test } from '@substrate-system/tapzero'
import type { Test } from '@substrate-system/tapzero'
import {
    gzipSync,
    zlibSync,
    deflateSync,
    Gunzip,
    Inflate,
    Unzlib,
    Decompress,
    strToU8,
    strFromU8
} from 'fflate'

// Regression coverage for TODO/callback.md.
//
// 0.8.3's sync-flush work made the streaming decompressors emit
// zero-length chunks with `final: false`. Those carry no data and no
// end-of-stream signal, so a consumer that keeps only the most recent
// chunk -- which is what @arethetypeswrong/core did -- ends up with an
// empty buffer. The rule enforced here is: a zero-length chunk is only
// ever emitted together with `final: true`.

type Emission = { len:number; final:boolean }

const record = (chunks:Uint8Array[], log:Emission[]) => {
    return (c:Uint8Array, final:boolean) => {
        chunks.push(c)
        log.push({ len:c.length, final })
    }
}

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

// The three invariants every streaming decompressor must hold, checked
// against one recorded run.
const assertContract = (t:Test, log:Emission[], label:string) => {
    const empties = log.filter(e => e.len === 0 && !e.final)
    t.equal(
        empties.length,
        0,
        label + ': no zero-length chunk is emitted with final false'
    )
    t.equal(
        log.filter(e => e.final).length,
        1,
        label + ': exactly one callback fires with final true'
    )
    t.equal(
        log[log.length - 1].final,
        true,
        label + ': the final callback is the last one'
    )
}

const HELLO = 'hello world'

test('the callback.md repro: Gunzip over a small gzip payload', t => {
    const gz = gzipSync(strToU8(HELLO))
    const chunks:Uint8Array[] = []
    const log:Emission[] = []
    new Gunzip(record(chunks, log)).push(gz, true)

    assertContract(t, log, 'Gunzip')
    t.equal(strFromU8(join(chunks)), HELLO, 'Gunzip: data round-trips')
})

test('no empty non-final chunks from Inflate, Unzlib or Decompress', t => {
    const raw = strToU8(HELLO)
    const cases:[string, () => { chunks:Uint8Array[]; log:Emission[] }][] = [
        ['Inflate', () => {
            const chunks:Uint8Array[] = []
            const log:Emission[] = []
            new Inflate(record(chunks, log)).push(deflateSync(raw), true)
            return { chunks, log }
        }],
        ['Unzlib', () => {
            const chunks:Uint8Array[] = []
            const log:Emission[] = []
            new Unzlib(record(chunks, log)).push(zlibSync(raw), true)
            return { chunks, log }
        }],
        ['Decompress/gzip', () => {
            const chunks:Uint8Array[] = []
            const log:Emission[] = []
            new Decompress(record(chunks, log)).push(gzipSync(raw), true)
            return { chunks, log }
        }],
        ['Decompress/zlib', () => {
            const chunks:Uint8Array[] = []
            const log:Emission[] = []
            new Decompress(record(chunks, log)).push(zlibSync(raw), true)
            return { chunks, log }
        }],
        ['Decompress/raw', () => {
            const chunks:Uint8Array[] = []
            const log:Emission[] = []
            new Decompress(record(chunks, log)).push(deflateSync(raw), true)
            return { chunks, log }
        }]
    ]

    for (const [label, run] of cases) {
        const { chunks, log } = run()
        assertContract(t, log, label)
        t.equal(strFromU8(join(chunks)), HELLO, label + ': data round-trips')
    }
})

test('multi-member gzip still emits one final chunk and no empties', t => {
    const a = gzipSync(strToU8('first member '))
    const b = gzipSync(strToU8('second member'))
    const both = new Uint8Array(a.length + b.length)
    both.set(a)
    both.set(b, a.length)

    const chunks:Uint8Array[] = []
    const log:Emission[] = []
    const members:number[] = []
    const gu = new Gunzip(record(chunks, log))
    gu.onmember = off => members.push(off)
    gu.push(both, true)

    assertContract(t, log, 'Gunzip multi-member')
    t.equal(
        strFromU8(join(chunks)),
        'first member second member',
        'Gunzip multi-member: both members are decompressed'
    )
    t.equal(members.length, 1, 'the second member is still reported')
})

test('an empty payload still gets exactly one final callback', t => {
    const chunks:Uint8Array[] = []
    const log:Emission[] = []
    new Gunzip(record(chunks, log)).push(gzipSync(new Uint8Array(0)), true)

    assertContract(t, log, 'Gunzip empty payload')
    t.equal(log.length, 1, 'end of stream is not lost when there is no data')
    t.equal(join(chunks).length, 0, 'and no bytes are invented')
})

test('chunked pushes still emit real bytes before the final chunk', t => {
    // Suppression must be limited to EMPTY chunks: intermediate chunks
    // that carry bytes -- which is what a sync flush boundary produces --
    // still have to reach the consumer as they arrive, not be held back
    // until the end. A payload well over the 32KiB window guarantees the
    // inflater has output to hand back before the last push.
    const src = new Uint8Array(200000)
    for (let i = 0; i < src.length; ++i) src[i] = (i * 7 + (i >> 5)) & 0xff
    const gz = gzipSync(src)

    const chunks:Uint8Array[] = []
    const log:Emission[] = []
    const gu = new Gunzip(record(chunks, log))
    const step = 4096
    for (let i = 0; i < gz.length; i += step) {
        const end = Math.min(i + step, gz.length)
        gu.push(gz.subarray(i, end), end === gz.length)
    }

    assertContract(t, log, 'Gunzip chunked')
    const midway = log.slice(0, -1).filter(e => e.len > 0)
    t.ok(
        midway.length > 0,
        'non-empty intermediate chunks are still delivered as they arrive'
    )

    const out = join(chunks)
    t.equal(out.length, src.length, 'chunked output has the right length')
    let same = true
    for (let i = 0; i < src.length; ++i) {
        if (out[i] !== src[i]) { same = false; break }
    }
    t.ok(same, 'chunked output matches the source bytes')
})
