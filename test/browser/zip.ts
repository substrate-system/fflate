import { test } from '@substrate-system/tapzero'
import type { Fflate } from './util.js'
import {
    fixtures,
    eq,
    concat,
    largeCompressible,
    largeIncompressible,
    withTimeout,
    firstDiff
} from './util.js'

// fflate stamps wall-clock time into every ZIP header
// (src/index.ts:2917) and the DOS time field has 2-second resolution,
// so two archives built moments apart differ at the mod-time offsets
// whenever the calls straddle a tick. Pin it on both sides of every
// byte-identity comparison.
const MTIME = new Date('2020-01-01T00:00:00Z')

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

// A filename in cp866 (the code page Russian Windows tools write),
// the ASCII placeholder that reserves exactly the same number of bytes
// for it, and two control names.
const CP866_NAME = 'привет.txt'
const CP866_PLACEHOLDER = 'PRIVET.TXT'
const UTF8_NAME = 'документ.txt'
const ASCII_NAME = 'plain.txt'

// cp866 maps every byte to exactly one character, so the decoder can
// be inverted to get an encoder. Building it this way keeps the
// expected bytes honest: they come from the same table the assertions
// decode with, rather than from a hand-copied literal.
function cp866Encode (s:string):Uint8Array {
    const dec = new TextDecoder('cp866')
    const rev = new Map<string, number>()
    for (let i = 0; i < 256; ++i) {
        rev.set(dec.decode(new Uint8Array([i])), i)
    }
    const out = new Uint8Array(s.length)
    for (let i = 0; i < s.length; ++i) {
        const b = rev.get(s[i])
        if (b == null) throw new Error('not in cp866: ' + s[i])
        out[i] = b
    }
    return out
}

// How fflate decodes a name whose UTF-8 bit is clear: byte-for-byte
// latin1, which round-trips losslessly but reads as mojibake.
function latin1 (bytes:Uint8Array):string {
    let s = ''
    for (let i = 0; i < bytes.length; ++i) s += String.fromCharCode(bytes[i])
    return s
}

// An archive with one entry whose name is cp866 and whose UTF-8 flag
// is clear. fflate's own zipSync always writes UTF-8 and sets bit 11
// for any non-ASCII name, so it cannot produce this directly: the
// entry goes in under an ASCII placeholder of the same byte length --
// which leaves bit 11 clear -- and the name bytes are then patched in
// place in both the local and the central header. Lengths are
// identical, so nothing else in the archive moves.
function cp866Fixture (f:Fflate):Uint8Array {
    const content = new TextEncoder().encode('hello world')
    const zipped = f.zipSync({
        [CP866_PLACEHOLDER]:content,
        [UTF8_NAME]:content,
        [ASCII_NAME]:content
    }, { level:0, mtime:MTIME })
    const to = cp866Encode(CP866_NAME)
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
        if (name === CP866_PLACEHOLDER) {
            if (fnl !== to.length) throw new Error('length mismatch')
            const lo = v.getUint32(o + 42, true)
            zipped.set(to, o + 46)
            zipped.set(to, lo + 30)
        }
        o += 46 + fnl + efl + cml
    }
    return zipped
}

// opts is required, not optional: zip() only substitutes {} when the
// third argument is absent (src/index.ts:3399), so passing an explicit
// undefined reaches fltn() and throws on `p.comment`.
function zipAsync (
    f:Fflate,
    data:Record<string, Uint8Array>,
    opts:Record<string, unknown>
):Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
        f.zip(data, opts, (err, result) => {
            if (err) reject(err)
            else resolve(result)
        })
    })
}

function unzipAsync (
    f:Fflate,
    data:Uint8Array
):Promise<Record<string, Uint8Array>> {
    return new Promise((resolve, reject) => {
        f.unzip(data, (err, result) => {
            if (err) reject(err)
            else resolve(result)
        })
    })
}

export function zipSuite (f:Fflate, label:string) {
    // Fixtures are regenerated per access, not shared. The async ZIP
    // APIs transfer their input buffers to a worker, which detaches
    // them -- a shared array would be emptied for every later test.
    const fx = () => fixtures()

    // Test basic zip/unzip round-trip
    test(`${label} > zipSync/unzipSync multi-file`, async t => {
        const compressible = fx()[0]
        const random = fx()[1]
        const text = fx()[2]
        const empty = fx()[3]

        const archive = f.zipSync({
            'file1.bin':compressible,
            'subdir/file2.bin':random,
            'text.txt':text,
            'empty.txt':empty
        })

        const unzipped = f.unzipSync(archive)

        t.ok(eq(unzipped['file1.bin'], compressible),
            'file1.bin content preserved')
        t.ok(eq(unzipped['subdir/file2.bin'], random),
            'nested file2.bin content preserved')
        t.ok(eq(unzipped['text.txt'], text),
            'text.txt content preserved')
        t.ok(eq(unzipped['empty.txt'], empty),
            'empty.txt content preserved')
    })

    // Test async zip/unzip against sync with small fixtures first
    test(`${label} > zip async vs sync small`, async t => {
        const compressible = fx()[0]
        const random = fx()[1]

        const syncArchive = f.zipSync({
            'file1.bin':compressible,
            'file2.bin':random
        }, { mtime:MTIME })

        try {
            const asyncArchive = await withTimeout(
                () => zipAsync(f, {
                    'file1.bin':compressible,
                    'file2.bin':random
                }, { mtime:MTIME })
            )
            const diff = firstDiff(asyncArchive, syncArchive)
            t.ok(diff === null, diff ?
        `async zip differs: ${diff}` :
                'async zip equals sync zip')
        } catch (e) {
            t.ok(false, `async zip failed: ${e}`)
        }
    })

    // Test async zip with large fixture to ensure worker is spawned
    test(`${label} > zip async vs sync large`, async t => {
        const largeData = largeCompressible()

        // Guard the premise. zip() compresses on the main thread when
        // originalSize < 160000 (src/index.ts:3466), and the archives
        // still match byte for byte when it does -- so without this a
        // shrunken fixture would silently remove the zip arm of the
        // wcln gate while the test stayed green.
        t.ok(largeData.length >= 160000,
            'zip worker branch needs originalSize >= 160000, ' +
      `got ${largeData.length}`)

        const syncArchive = f.zipSync({
            'large.bin':largeData
        }, { mtime:MTIME })

        try {
            const asyncArchive = await withTimeout(
                () => zipAsync(f, {
                    'large.bin':largeData
                }, { mtime:MTIME })
            )
            const diff = firstDiff(asyncArchive, syncArchive)
            t.ok(diff === null, diff ?
        `async zip large differs: ${diff}` :
                'async zip large equals sync zip')
        } catch (e) {
            t.ok(false, `async zip large failed: ${e}`)
        }
    })

    // Test async unzip with small fixtures first
    test(`${label} > unzip async vs sync small`, async t => {
        const compressible = fx()[0]
        const random = fx()[1]

        const archive = f.zipSync({
            'file1.bin':compressible,
            'file2.bin':random
        })

        const syncUnzip = f.unzipSync(archive)

        try {
            const asyncUnzip = await withTimeout(
                () => unzipAsync(f, archive)
            )

            for (const key of Object.keys(syncUnzip)) {
                t.ok(key in asyncUnzip, `key ${key} in async result`)
                if (key in asyncUnzip) {
                    t.ok(eq(asyncUnzip[key], syncUnzip[key]),
            `async unzip key ${key} matches sync`)
                }
            }
        } catch (e) {
            t.ok(false, `async unzip failed: ${e}`)
        }
    })

    // Test async unzip with a large COMPRESSIBLE fixture, so that
    // unzip() takes its worker branch. It does so only when
    // originalSize >= 524288 AND compressedSize <= 0.8 * originalSize,
    // which an incompressible fixture would fail on the second clause.
    test(`${label} > unzip async vs sync large`, async t => {
        const largeData = largeCompressible()

        const archive = f.zipSync({
            'large.bin':largeData
        })

        // Guard the premise. zipSync and deflateSync use the same
        // defaults, so this is the compressed size of the entry itself.
        // Without these the test still passes on a fixture that never
        // leaves the main thread.
        const su = largeData.length
        const sc = f.deflateSync(largeData).length
        t.ok(su >= 524288,
      `unzip worker branch needs originalSize >= 524288, got ${su}`)
        t.ok(sc <= 0.8 * su,
            'unzip worker branch needs compressed <= 0.8 * original, ' +
      `got ${sc} vs ${0.8 * su}`)

        const syncUnzip = f.unzipSync(archive)

        try {
            const asyncUnzip = await withTimeout(
                () => unzipAsync(f, archive)
            )

            for (const key of Object.keys(syncUnzip)) {
                t.ok(key in asyncUnzip, `key ${key} in async result`)
                if (key in asyncUnzip) {
                    t.ok(eq(asyncUnzip[key], syncUnzip[key]),
            `async unzip key ${key} matches sync`)
                }
            }
        } catch (e) {
            t.ok(false, `async unzip large failed: ${e}`)
        }
    })

    // Test streaming Zip with ZipDeflate
    test(
    `${label} > Zip streaming ZipDeflate`,
    async t => {
        const compressible = fx()[0]
        const random = fx()[1]

        const zipped:Uint8Array[] = []
        const zipper = new f.Zip((err, chunk, _final) => {
            if (!err) {
                zipped.push(new Uint8Array(chunk))
            }
        })

        // add() installs the entry's ondata handler, so it must come
        // before any push -- otherwise ZipPassThrough.push errors with
        // 'no stream handler'.
        const file1 = new f.ZipDeflate('file1.bin')
        zipper.add(file1)
        file1.push(compressible.slice(), true)

        const file2 = new f.ZipDeflate('file2.bin')
        zipper.add(file2)
        file2.push(random.slice(), true)

        zipper.end()

        const unzipped = f.unzipSync(concat(zipped))
        t.ok(eq(unzipped['file1.bin'], compressible),
            'streaming Zip file1 preserved')
        t.ok(eq(unzipped['file2.bin'], random),
            'streaming Zip file2 preserved')
    }
    )

    // Test ZipPassThrough (stored, not deflated)
    test(
    `${label} > ZipPassThrough stored entry`,
    async t => {
        const compressible = fx()[0]

        const zipped:Uint8Array[] = []
        const zipper = new f.Zip((err, chunk, _final) => {
            if (!err) {
                zipped.push(new Uint8Array(chunk))
            }
        })

        // Add as pass-through (stored)
        const file = new f.ZipPassThrough('stored.bin')
        zipper.add(file)
        file.push(compressible.slice(), true)

        zipper.end()

        const unzipped = f.unzipSync(concat(zipped))
        t.ok(eq(unzipped['stored.bin'], compressible),
            'ZipPassThrough stored entry preserved')
    }
    )

    // Test AsyncZipDeflate streaming
    test(
    `${label} > Zip streaming AsyncZipDeflate`,
    async t => {
        const compressible = fx()[0]

        const zipped:Uint8Array[] = []
        let zipError:Error|null = null
        let zipSettled = false
        // Zip delivers worker errors with `final` false
        // (src/index.ts:3311), so settle on either. Settling only on
        // `final` would spin to the withTimeout deadline on an error
        // instead of failing with its text.
        const zipper = new f.Zip((err, chunk, final) => {
            if (err) {
                zipError = err
                zipSettled = true
                return
            }
            zipped.push(new Uint8Array(chunk))
            if (final) zipSettled = true
        })

        // Push a COPY: the async entry transfers its input buffer to
        // the worker, which detaches it. Pushing the shared fixture
        // directly would empty it for every later test in this file.
        const file = new f.AsyncZipDeflate('async.bin')
        zipper.add(file)
        file.push(compressible.slice(), true)

        zipper.end()

        try {
            await withTimeout(() =>
                new Promise<void>(resolve => {
                    const check = () => {
                        if (zipSettled) {
                            resolve()
                        } else {
                            setTimeout(check, 10)
                        }
                    }
                    check()
                })
            )

            t.ok(zipError === null, zipError ?
          `AsyncZipDeflate error: ${zipError}` :
                'AsyncZipDeflate reported no error')

            const unzipped = f.unzipSync(concat(zipped))
            t.ok(eq(unzipped['async.bin'], compressible),
                'AsyncZipDeflate streaming preserved')
        } catch (e) {
            t.ok(false, `AsyncZipDeflate failed: ${e}`)
        }
    }
    )

    // Test streaming Unzip via Unzip class
    test(
    `${label} > Unzip streaming basic`,
    async t => {
        const compressible = fx()[0]
        const random = fx()[1]

        const archive = f.zipSync({
            'file1.bin':compressible,
            'file2.bin':random
        })

        const entries:Record<string, Uint8Array[]> = {}

        const unzipper = new f.Unzip((file) => {
            if (!entries[file.name]) {
                entries[file.name] = []
            }
            file.ondata = (err, data, _final) => {
                if (!err) {
                    entries[file.name].push(new Uint8Array(data))
                }
            }
            // Nothing flows until start() is called.
            file.start()
        })

        // Entries written by zipSync are deflated, so Unzip needs a
        // decoder registered for them. Without this it only handles
        // stored entries and deflated ones yield no data at all.
        unzipper.register(f.UnzipInflate)

        // Push archive in chunks
        const chunkSize = Math.ceil(archive.length / 3)
        for (let i = 0; i < archive.length; i += chunkSize) {
            const end = Math.min(i + chunkSize, archive.length)
            const isLast = end === archive.length
            unzipper.push(archive.slice(i, end), isLast)
        }

        const reconstructed:Record<string, Uint8Array> = {}
        for (const [name, parts] of Object.entries(entries)) {
            reconstructed[name] = concat(parts)
        }

        t.ok(eq(reconstructed['file1.bin'], compressible),
            'streaming Unzip file1 preserved')
        t.ok(eq(reconstructed['file2.bin'], random),
            'streaming Unzip file2 preserved')
    }
    )

    // Streaming Unzip reads LOCAL file headers, which carry no external
    // file attributes, so isDirectory there can only come from the
    // trailing-slash convention. The archive below deliberately mixes a
    // directory entry, an empty file, and a normal file so that the
    // empty file is not misreported as a directory.
    test(
    `${label} > Unzip streaming isDirectory`,
    async t => {
        const archive = f.zipSync({
            'dir/':new Uint8Array(0),
            'dir/empty.bin':new Uint8Array(0),
            'dir/file.bin':fixtures()[0]
        }, { mtime:MTIME })

        const seen:Record<string, boolean> = {}
        const order:string[] = []

        const unzipper = new f.Unzip((file) => {
            order.push(file.name)
            seen[file.name] = file.isDirectory
            file.ondata = () => {}
            file.start()
        })
        unzipper.register(f.UnzipInflate)

        unzipper.push(archive, true)

        t.equal(order.length, 3, 'streaming Unzip saw every entry')
        t.equal(seen['dir/'], true, 'directory entry is a directory')
        t.equal(seen['dir/empty.bin'], false,
            'empty file is not a directory')
        t.equal(seen['dir/file.bin'], false,
            'normal file is not a directory')
    }
    )

    // ZIP has no encoding field for filenames. Bit 11 of the general
    // purpose bit flag is the only signal, and it only ever says
    // "UTF-8"; when it is clear, the bytes are in whatever code page
    // the writing tool used. Streaming extraction reads LOCAL file
    // headers, which carry the same flag and the same name bytes as
    // the central directory, so the raw bytes are available here too.
    test(
    `${label} > Unzip streaming rawName and nameIsUTF8`,
    async t => {
        const archive = cp866Fixture(f)

        const seen:Record<string, {
            rawName:Uint8Array;
            nameIsUTF8:boolean;
        }> = {}

        const unzipper = new f.Unzip((file) => {
            seen[file.name] = {
                // the bytes must survive the callback; copy so a later
                // push cannot change what was captured
                rawName:new Uint8Array(file.rawName),
                nameIsUTF8:file.nameIsUTF8
            }
            file.ondata = () => {}
            file.start()
        })
        unzipper.register(f.UnzipInflate)

        // Chunked, so the header lands across a push boundary at least
        // some of the time: the raw bytes come out of Unzip's internal
        // buffer, not out of the chunk the caller handed in.
        const chunkSize = Math.ceil(archive.length / 4)
        for (let i = 0; i < archive.length; i += chunkSize) {
            const end = Math.min(i + chunkSize, archive.length)
            unzipper.push(archive.slice(i, end), end === archive.length)
        }

        const mojibake = latin1(cp866Encode(CP866_NAME))
        const cp866 = seen[mojibake]
        t.ok(cp866, 'streaming Unzip saw the cp866 entry')
        t.equal(cp866.nameIsUTF8, false, 'the cp866 entry is not UTF-8')
        t.deepEqual(
            Array.from(cp866.rawName),
            Array.from(cp866Encode(CP866_NAME)),
            'rawName is the bytes exactly as stored'
        )
        t.equal(
            new TextDecoder('cp866').decode(cp866.rawName),
            CP866_NAME,
            'decoding rawName as cp866 yields the real filename'
        )

        t.equal(seen[UTF8_NAME].nameIsUTF8, true,
            'a non-ASCII name written by fflate is marked UTF-8')
        t.deepEqual(
            Array.from(seen[UTF8_NAME].rawName),
            Array.from(new TextEncoder().encode(UTF8_NAME)),
            'a UTF-8 name re-encodes to exactly the stored bytes')

        t.equal(seen[ASCII_NAME].nameIsUTF8, false,
            'an ASCII name has bit 11 clear, and needs no decoding')
        t.deepEqual(
            Array.from(seen[ASCII_NAME].rawName),
            Array.from(new TextEncoder().encode(ASCII_NAME)),
            'rawName matches name for an ASCII entry')
    }
    )

    // Test UnzipInflate: sync decompression during Unzip
    test(
    `${label} > UnzipInflate sync decompression`,
    async t => {
        const compressible = fx()[0]

        const archive = f.zipSync({
            'deflated.bin':compressible
        })

        const recovered:Record<string, Uint8Array> = {}
        // Poll on decodeSettled, which the error branch sets too, not on
        // recovered[name], which it never sets. Polling on the payload
        // would leave the assertion below stranded behind an await that
        // only resolves on success -- unreachable in the one case it
        // exists for, so the decoder error would surface as a bare
        // `timeout` with its text lost.
        let decodeErr:Error|null = null
        let decodeSettled = false

        // Drive the decoder through the public API: register it and
        // let Unzip construct it. Manually constructing the decoder
        // and calling file.start() conflicts -- start() instantiates
        // the REGISTERED constructor and throws without one.
        const unzipper = new f.Unzip((file) => {
            if (file.name === 'deflated.bin') {
                const parts:Uint8Array[] = []
                file.ondata = (err, data, final) => {
                    if (err) {
                        decodeErr = err
                        decodeSettled = true
                        return
                    }
                    parts.push(new Uint8Array(data))
                    if (final) {
                        recovered[file.name] = concat(parts)
                        decodeSettled = true
                    }
                }
                file.start()
            }
        })
        unzipper.register(f.UnzipInflate)

        const chunkSize = Math.ceil(archive.length / 3)
        for (let i = 0; i < archive.length; i += chunkSize) {
            const end = Math.min(i + chunkSize, archive.length)
            const isLast = end === archive.length
            unzipper.push(archive.slice(i, end), isLast)
        }

        try {
            await withTimeout(() =>
                new Promise<void>(resolve => {
                    const check = () => {
                        if (decodeSettled) {
                            resolve()
                        } else {
                            setTimeout(check, 10)
                        }
                    }
                    check()
                })
            )

            t.ok(decodeErr === null, decodeErr ?
          `decoder error: ${decodeErr}` :
                'decoder reported no error')

            const entry = recovered['deflated.bin']
            t.ok(eq(entry, compressible),
                'UnzipInflate sync decompression preserved')
        } catch (e) {
            t.ok(false, `UnzipInflate sync decompression failed: ${e}`)
        }
    }
    )

    // Test AsyncUnzipInflate: async decompression during Unzip
    // Uses small fixture first
    test(
    `${label} > AsyncUnzipInflate async decompression small`,
    async t => {
        const compressible = fx()[0]

        const archive = f.zipSync({
            'async-deflated.bin':compressible
        })

        const recovered:Record<string, Uint8Array> = {}
        // Poll on decodeSettled, which the error branch sets too, not on
        // recovered[name], which it never sets. Polling on the payload
        // would leave the assertion below stranded behind an await that
        // only resolves on success -- unreachable in the one case it
        // exists for, so the decoder error would surface as a bare
        // `timeout` with its text lost.
        let decodeErr:Error|null = null
        let decodeSettled = false

        // Drive the decoder through the public API: register it and
        // let Unzip construct it. Manually constructing the decoder
        // and calling file.start() conflicts -- start() instantiates
        // the REGISTERED constructor and throws without one.
        const unzipper = new f.Unzip((file) => {
            if (file.name === 'async-deflated.bin') {
                const parts:Uint8Array[] = []
                file.ondata = (err, data, final) => {
                    if (err) {
                        decodeErr = err
                        decodeSettled = true
                        return
                    }
                    parts.push(new Uint8Array(data))
                    if (final) {
                        recovered[file.name] = concat(parts)
                        decodeSettled = true
                    }
                }
                file.start()
            }
        })
        unzipper.register(f.AsyncUnzipInflate)

        const chunkSize = Math.ceil(archive.length / 3)
        for (let i = 0; i < archive.length; i += chunkSize) {
            const end = Math.min(i + chunkSize, archive.length)
            const isLast = end === archive.length
            unzipper.push(archive.slice(i, end), isLast)
        }

        try {
            await withTimeout(() =>
                new Promise<void>(resolve => {
                    const check = () => {
                        if (decodeSettled) {
                            resolve()
                        } else {
                            setTimeout(check, 10)
                        }
                    }
                    check()
                })
            )

            t.ok(decodeErr === null, decodeErr ?
          `decoder error: ${decodeErr}` :
                'decoder reported no error')

            const entry = recovered['async-deflated.bin']
            t.ok(
                eq(entry, compressible),
                'AsyncUnzipInflate async decompression preserved'
            )
        } catch (e) {
            t.ok(false, `AsyncUnzipInflate async decompression small failed: ${e}`)
        }
    }
    )

    // Test AsyncUnzipInflate with large incompressible fixture to
    // ensure async decoder actually spawns a worker. The compressed
    // size must exceed 320000 bytes for the async path to trigger.
    test(
    `${label} > AsyncUnzipInflate async decompression large`,
    async t => {
        const largeData = largeIncompressible()

        const archive = f.zipSync({
            'async-large.bin':largeData
        })

        const recovered:Record<string, Uint8Array> = {}
        // Poll on decodeSettled, which the error branch sets too, not on
        // recovered[name], which it never sets. Polling on the payload
        // would leave the assertion below stranded behind an await that
        // only resolves on success -- unreachable in the one case it
        // exists for, so the decoder error would surface as a bare
        // `timeout` with its text lost.
        let decodeErr:Error|null = null
        let decodeSettled = false
        let entrySize = -1

        // Unzip constructs the registered decoder itself and keeps the
        // instance private, so wrap AsyncUnzipInflate to observe which
        // branch it took. `terminate` is assigned only on the async
        // branch, so it is the one externally visible difference
        // between AsyncInflate and the synchronous Inflate fallback.
        //
        // A flag set inside the onfile callback cannot serve here: that
        // callback fires for the sync decoder too, so the assertion
        // could never fail.
        let sawAsyncBranch = false

        type Decoder = InstanceType<typeof f.AsyncUnzipInflate>

        class SpyUnzipInflate {
            static compression = 8
            private i:Decoder
            ondata:Decoder['ondata']
            terminate:Decoder['terminate']

            constructor (fn:string, sz?:number) {
                this.i = new f.AsyncUnzipInflate(fn, sz)
                this.i.ondata = (err, dat, final) => {
                    this.ondata(err, dat, final)
                }
                if (this.i.terminate) {
                    sawAsyncBranch = true
                    this.terminate = this.i.terminate
                }
            }

            push (chunk:Uint8Array, final:boolean) {
                this.i.push(chunk, final)
            }
        }

        const unzipper = new f.Unzip((file) => {
            if (file.name === 'async-large.bin') {
                entrySize = file.size
                const parts:Uint8Array[] = []
                file.ondata = (err, data, final) => {
                    if (err) {
                        decodeErr = err
                        decodeSettled = true
                        return
                    }
                    parts.push(new Uint8Array(data))
                    if (final) {
                        recovered[file.name] = concat(parts)
                        decodeSettled = true
                    }
                }
                file.start()
            }
        })
        unzipper.register(SpyUnzipInflate)

        const chunkSize = Math.ceil(archive.length / 3)
        for (let i = 0; i < archive.length; i += chunkSize) {
            const end = Math.min(i + chunkSize, archive.length)
            const isLast = end === archive.length
            unzipper.push(archive.slice(i, end), isLast)
        }

        try {
            await withTimeout(() =>
                new Promise<void>(resolve => {
                    const check = () => {
                        if (decodeSettled) {
                            resolve()
                        } else {
                            setTimeout(check, 10)
                        }
                    }
                    check()
                })
            )

            // Unzip passes the COMPRESSED size as `sz`, so an
            // over-compressible fixture keeps this under 320000 and the
            // async branch is never reached.
            t.ok(entrySize >= 320000,
                'AsyncUnzipInflate async branch needs compressed size ' +
          `>= 320000, got ${entrySize}`)
            t.ok(sawAsyncBranch,
                'AsyncUnzipInflate large took the AsyncInflate branch')

            t.ok(decodeErr === null, decodeErr ?
          `decoder error: ${decodeErr}` :
                'decoder reported no error')

            const entry = recovered['async-large.bin']
            t.ok(
                eq(entry, largeData),
                'AsyncUnzipInflate large async decompression preserved'
            )
        } catch (e) {
            t.ok(false, `AsyncUnzipInflate async decompression large failed: ${e}`)
        }
    }
    )

    // The other half of AC3.4's "each is readable by the other's
    // decoder". The tests above read a sync-produced archive with the
    // async decoder; this reads an async-produced one with unzipSync.
    test(`${label} > unzipSync reads async-produced archive`, async t => {
        const large = largeCompressible()
        const small = fx()[1]

        // At least one entry must clear zip()'s 160000 gate
        // (src/index.ts:3466), or the archive is built entirely on the
        // main thread and this proves nothing about the worker path --
        // which is the only hazard this file exists to catch.
        t.ok(large.length >= 160000,
            'zip worker branch needs originalSize >= 160000, ' +
      `got ${large.length}`)

        try {
            const asyncArchive = await withTimeout(
                () => zipAsync(f, {
                    'async-large.bin':large,
                    'async-small.bin':small
                }, { mtime:MTIME })
            )

            // Read with sync decoder
            const syncUnzip = f.unzipSync(asyncArchive)

            t.ok(eq(syncUnzip['async-large.bin'], large),
                'unzipSync reads async-produced worker entry')
            t.ok(eq(syncUnzip['async-small.bin'], small),
                'unzipSync reads async-produced main-thread entry')
        } catch (e) {
            t.ok(false, `unzipSync reads async-produced archive failed: ${e}`)
        }
    })

    // One archive carrying all three entry types at once. The
    // per-type tests above each build a single-file archive, which
    // never exercises stored and deflated entries sharing a central
    // directory, nor an async entry finishing out of order.
    test(
    `${label} > Zip streaming mixed entry types`,
    async t => {
        const compressible = fx()[0]
        const random = fx()[1]
        const text = fx()[2]

        const zipped:Uint8Array[] = []
        let zipError:Error|null = null
        let zipSettled = false
        // Settle on either arm -- see Zip streaming AsyncZipDeflate.
        const zipper = new f.Zip((err, chunk, final) => {
            if (err) {
                zipError = err
                zipSettled = true
                return
            }
            zipped.push(new Uint8Array(chunk))
            if (final) zipSettled = true
        })

        // Mix different entry types
        const deflateEntry = new f.ZipDeflate('deflated.bin')
        zipper.add(deflateEntry)
        deflateEntry.push(compressible.slice(), true)

        const storedEntry = new f.ZipPassThrough('stored.txt')
        zipper.add(storedEntry)
        storedEntry.push(text.slice(), true)

        const asyncEntry = new f.AsyncZipDeflate('async.bin')
        zipper.add(asyncEntry)
        asyncEntry.push(random.slice(), true)

        zipper.end()

        try {
            await withTimeout(() =>
                new Promise<void>(resolve => {
                    const check = () => {
                        if (zipSettled) {
                            resolve()
                        } else {
                            setTimeout(check, 10)
                        }
                    }
                    check()
                })
            )

            t.ok(zipError === null, zipError ?
          `mixed archive error: ${zipError}` :
                'mixed archive reported no error')

            const unzipped = f.unzipSync(concat(zipped))

            t.ok(eq(unzipped['deflated.bin'], compressible),
                'mixed archive deflated.bin preserved')
            t.ok(eq(unzipped['stored.txt'], text),
                'mixed archive stored.txt preserved')
            t.ok(eq(unzipped['async.bin'], random),
                'mixed archive async.bin preserved')
        } catch (e) {
            t.ok(false, `Zip mixed entry types failed: ${e}`)
        }
    }
    )

    // The DOS mod-time word lives in the central directory header, which
    // is the header zh() parses for unzipSync. Two entries with
    // DIFFERENT times is the point: a single shared read would pass an
    // equal-times test.
    test(`${label} > unzipSync filter receives per-entry mtime`, async t => {
        const compressible = fx()[0]
        const random = fx()[1]

        const archive = f.zipSync({
            'a.bin':[compressible, { mtime:MTIME_A }],
            'b.bin':[random, { mtime:MTIME_B }]
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
            'min.bin':[compressible, { mtime:MTIME_MIN }],
            'max.bin':[random, { mtime:MTIME_MAX }]
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

    // Same central-header path as unzipSync, but through the async
    // control flow. The unzipAsync() helper at the top of this file uses
    // the two-argument form, which cannot carry a filter, so this test
    // drives f.unzip directly.
    test(`${label} > unzip filter receives per-entry mtime`, async t => {
        const compressible = fx()[0]
        const random = fx()[1]

        const archive = f.zipSync({
            'a.bin':[compressible, { mtime:MTIME_A }],
            'b.bin':[random, { mtime:MTIME_B }]
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

    // Unzip parses LOCAL file headers, not the central directory, so
    // this exercises a different offset (+10) and a different code path
    // from the two filter tests above.
    test(`${label} > Unzip stream exposes per-entry mtime`, async t => {
        const compressible = fx()[0]
        const random = fx()[1]

        const archive = f.zipSync({
            'a.bin':[compressible, { mtime:MTIME_A }],
            'b.bin':[random, { mtime:MTIME_B }]
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
}
