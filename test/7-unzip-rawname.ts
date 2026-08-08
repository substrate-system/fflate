import { test } from '@substrate-system/tapzero'
import { zipSync, unzipSync, unzip, type UnzipFileInfo } from 'fflate'

// ZIP stores filenames as bytes with no encoding field. Bit 11 of the
// general purpose bit flag is the only signal, and it only ever says
// "UTF-8"; when it is clear the bytes are in whatever code page the
// writing tool used. fflate cannot guess that, so it exposes the raw
// bytes and lets the caller pick a TextDecoder.
//
// fflate's own zipSync always writes UTF-8 and sets bit 11 whenever the
// name is non-ASCII, so it cannot produce a cp866 fixture directly. The
// archive below is written with ASCII placeholder names of exactly the
// right byte length -- which leaves bit 11 clear -- and the filename
// bytes are then patched in place in both the local and the central
// header. The result is a real ZIP that any tool would accept, and is
// byte-for-byte what a Russian Windows tool would have produced.

// cp866 maps every byte to exactly one character, so the decoder can be
// inverted to get an encoder. Building it this way keeps the expected
// bytes honest: they come from the same table the assertion decodes
// with, rather than from a hand-copied literal.
const cp866Encode = (s:string) => {
    const dec = new TextDecoder('cp866')
    const rev = new Map<string, number>()
    for (let i = 0; i < 256; ++i) {
        rev.set(dec.decode(new Uint8Array([i])), i)
    }
    const out = new Uint8Array(s.length)
    for (let i = 0; i < s.length; ++i) {
        const b = rev.get(s[i])
        if (b == null) throw new Error('not representable in cp866: ' + s[i])
        out[i] = b
    }
    return out
}

// The name we want in the archive, and the ASCII placeholder that
// reserves exactly the same number of bytes for it.
const CP866_NAME = 'привет.txt'
const PLACEHOLDER = 'PRIVET.TXT'

const UTF8_NAME = 'документ.txt'
const ASCII_NAME = 'plain.txt'

const CONTENT = new TextEncoder().encode('hello world')

// Overwrite the filename bytes of one entry, in both the central
// directory record and the local file header it points at. Lengths are
// identical, so nothing else in the archive moves.
const patchName = (zipped:Uint8Array, from:string, to:Uint8Array) => {
    const v = new DataView(
        zipped.buffer,
        zipped.byteOffset,
        zipped.byteLength
    )
    let e = zipped.length - 22
    while (v.getUint32(e, true) !== 0x06054b50) --e
    const n = v.getUint16(e + 8, true)
    let o = v.getUint32(e + 16, true)
    const dec = new TextDecoder()
    for (let i = 0; i < n; ++i) {
        const fnl = v.getUint16(o + 28, true)
        const efl = v.getUint16(o + 30, true)
        const cml = v.getUint16(o + 32, true)
        const name = dec.decode(zipped.subarray(o + 46, o + 46 + fnl))
        if (name === from) {
            if (fnl !== to.length) throw new Error('length mismatch')
            const lo = v.getUint32(o + 42, true)
            zipped.set(to, o + 46)
            zipped.set(to, lo + 30)
        }
        o += 46 + fnl + efl + cml
    }
    return zipped
}

const buildFixture = () => {
    const zipped = zipSync({
        [PLACEHOLDER]:CONTENT,
        [UTF8_NAME]:CONTENT,
        [ASCII_NAME]:CONTENT
    }, { level:0 })
    return patchName(zipped, PLACEHOLDER, cp866Encode(CP866_NAME))
}

// The key unzipSync will use for the cp866 entry: fflate decodes a
// non-UTF-8 name as latin1, which round-trips the bytes but reads as
// mojibake.
const latin1Name = () => {
    const bytes = cp866Encode(CP866_NAME)
    let s = ''
    for (let i = 0; i < bytes.length; ++i) s += String.fromCharCode(bytes[i])
    return s
}

const collect = (zipped:Uint8Array) => {
    const seen:Record<string, UnzipFileInfo> = {}
    unzipSync(zipped, {
        filter:file => {
            seen[file.name] = {
                ...file,
                // the filter argument is reused across entries in
                // principle; copy the view so later entries cannot
                // change what we captured
                rawName:file.rawName.slice()
            }
            return true
        }
    })
    return seen
}

test('rawName carries the undecoded filename bytes', t => {
    const zipped = buildFixture()
    const seen = collect(zipped)
    const info = seen[latin1Name()]
    t.ok(info, 'the cp866 entry is present under its latin1 name')
    t.ok(info.rawName instanceof Uint8Array, 'rawName is a Uint8Array')
    t.deepEqual(
        Array.from(info.rawName),
        Array.from(cp866Encode(CP866_NAME)),
        'rawName is the bytes exactly as stored'
    )
    t.equal(
        new TextDecoder('cp866').decode(info.rawName),
        CP866_NAME,
        'decoding rawName as cp866 yields the real filename'
    )
})

test('nameIsUTF8 reflects bit 11 of the general purpose bit flag', t => {
    const seen = collect(buildFixture())
    t.equal(
        seen[latin1Name()].nameIsUTF8,
        false,
        'a cp866 entry is not marked UTF-8'
    )
    t.equal(
        seen[UTF8_NAME].nameIsUTF8,
        true,
        'a non-ASCII name written by fflate is marked UTF-8'
    )
    t.equal(
        seen[ASCII_NAME].nameIsUTF8,
        false,
        'an ASCII name has bit 11 clear, and needs no decoding anyway'
    )
})

test('name keeps its existing decoding behavior', t => {
    const seen = collect(buildFixture())
    t.equal(
        seen[UTF8_NAME].name,
        UTF8_NAME,
        'a UTF-8 name still decodes as UTF-8'
    )
    t.equal(
        seen[ASCII_NAME].name,
        ASCII_NAME,
        'an ASCII name is unchanged'
    )
    t.deepEqual(
        Array.from(seen[ASCII_NAME].rawName),
        Array.from(new TextEncoder().encode(ASCII_NAME)),
        'rawName matches name for an ASCII entry'
    )
    t.deepEqual(
        Array.from(new TextEncoder().encode(seen[UTF8_NAME].name)),
        Array.from(seen[UTF8_NAME].rawName),
        'a UTF-8 name re-encodes to exactly the stored bytes'
    )
})

test('unzip exposes the same fields as unzipSync', t => {
    return new Promise<void>((resolve, reject) => {
        const zipped = buildFixture()
        const seen:Record<string, UnzipFileInfo> = {}
        unzip(zipped, {
            filter:file => {
                seen[file.name] = { ...file, rawName:file.rawName.slice() }
                return true
            }
        }, err => {
            if (err) return reject(err)
            try {
                const info = seen[latin1Name()]
                t.equal(info.nameIsUTF8, false, 'async: not marked UTF-8')
                t.equal(
                    new TextDecoder('cp866').decode(info.rawName),
                    CP866_NAME,
                    'async: rawName decodes as cp866'
                )
                t.equal(
                    seen[UTF8_NAME].nameIsUTF8,
                    true,
                    'async: the UTF-8 entry is marked UTF-8'
                )
                resolve()
            } catch (e) {
                reject(e)
            }
        })
    })
})
