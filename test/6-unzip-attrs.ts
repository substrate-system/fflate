import { test } from '@substrate-system/tapzero'
import { zipSync, unzipSync, unzip, type UnzipFileInfo } from 'fflate'

// fflate's own zipSync leaves the external file attributes dword at 0,
// which is exactly the case the trailing-slash fallback exists for. To
// cover the DOS and unix attribute bits too we write a normal archive
// and then patch the central directory in place, so the fixture stays a
// real ZIP that any tool would accept.

const S_IFDIR = 0o040000
const S_IFREG = 0o100000
const DOS_DIR = 0x10

const patchAttrs = (zipped:Uint8Array, attrs:Record<string, number>) => {
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
        if (name in attrs) v.setUint32(o + 38, attrs[name], true)
        o += 46 + fnl + efl + cml
    }
    return zipped
}

const CONTENT = new TextEncoder().encode('hello world')

// `attrs` is read with an unsigned 32-bit read, so any unix st_mode
// shifted into the high half comes back above 2^31. Expectations have to
// be unsigned too: `0o100644 << 16` is negative in JS.
const mode = (m:number) => (m << 16) >>> 0

// name -> [data, external attributes]
const ENTRIES:Record<string, [Uint8Array, number]> = {
    // unix S_IFDIR only, no trailing slash
    unixdir:[new Uint8Array(0), mode(S_IFDIR | 0o755)],
    // MS-DOS directory bit only, no trailing slash
    dosdir:[new Uint8Array(0), DOS_DIR],
    // no attributes at all: only the trailing slash marks it
    'slashdir/':[new Uint8Array(0), 0],
    // an empty FILE -- the case that is indistinguishable from a
    // directory without this metadata
    'empty.txt':[new Uint8Array(0), mode(S_IFREG | 0o644)],
    // an ordinary file with content
    'file.txt':[CONTENT, mode(S_IFREG | 0o644)]
}

const buildFixture = () => {
    const raw:Record<string, Uint8Array> = {}
    const attrs:Record<string, number> = {}
    for (const k in ENTRIES) {
        raw[k] = ENTRIES[k][0]
        attrs[k] = ENTRIES[k][1]
    }
    return patchAttrs(zipSync(raw, { level:0 }), attrs)
}

const EXPECTED_DIR:Record<string, boolean> = {
    unixdir:true,
    dosdir:true,
    'slashdir/':true,
    'empty.txt':false,
    'file.txt':false
}

test('unzipSync filter receives attrs and isDirectory', t => {
    const zipped = buildFixture()
    const seen:Record<string, UnzipFileInfo> = {}
    unzipSync(zipped, {
        filter:file => {
            seen[file.name] = file
            return true
        }
    })

    t.equal(
        Object.keys(seen).length,
        Object.keys(ENTRIES).length,
        'the filter should see every entry in the archive'
    )

    for (const name in ENTRIES) {
        t.equal(
            seen[name].attrs,
            ENTRIES[name][1],
            `${name} should report its raw external file attributes`
        )
        t.equal(
            seen[name].isDirectory,
            EXPECTED_DIR[name],
            `${name} isDirectory should be ${EXPECTED_DIR[name]}`
        )
    }
})

test('isDirectory tells a directory apart from an empty file', t => {
    const zipped = buildFixture()
    const seen:Record<string, boolean> = {}
    const extracted = unzipSync(zipped, {
        filter:file => {
            seen[file.name] = file.isDirectory
            return true
        }
    })

    t.equal(
        extracted['empty.txt'].length,
        0,
        'empty.txt should decompress to zero bytes'
    )
    t.equal(
        extracted.unixdir.length,
        0,
        'the directory entry should also be zero bytes'
    )
    t.notEqual(
        seen['empty.txt'],
        seen.unixdir,
        'a zero-byte file and a zero-byte directory must not look alike'
    )
})

test('unzip filter receives attrs and isDirectory', t => {
    return new Promise<void>((resolve, reject) => {
        const zipped = buildFixture()
        const seen:Record<string, UnzipFileInfo> = {}
        unzip(zipped, {
            filter:file => {
                seen[file.name] = file
                return true
            }
        }, err => {
            if (err) return reject(err)
            try {
                for (const name in ENTRIES) {
                    t.equal(
                        seen[name].attrs,
                        ENTRIES[name][1],
                        `${name} attrs should match through unzip`
                    )
                    t.equal(
                        seen[name].isDirectory,
                        EXPECTED_DIR[name],
                        `${name} isDirectory should match through unzip`
                    )
                }
                resolve()
            } catch (e) {
                reject(e)
            }
        })
    })
})
