import { test } from '@substrate-system/tapzero'
import type { Fflate } from './util.js'
import { fixtures, eq, withTimeout, firstDiff } from './util.js'

// gzip headers carry a 1-second-resolution mtime, so any byte-identity
// comparison must pin it on both sides. See the gzip sync/async test.
const MTIME = new Date('2020-01-01T00:00:00Z')

function deflateAsync (
    f:Fflate,
    data:Uint8Array
):Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
        f.deflate(data, (err, result) => {
            if (err) reject(err)
            else {
                if (err !== null) reject(new Error('err not null'))
                else resolve(result!)
            }
        })
    })
}

function inflateAsync (
    f:Fflate,
    data:Uint8Array
):Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
        f.inflate(data, (err, result) => {
            if (err) reject(err)
            else {
                if (err !== null) reject(new Error('err not null'))
                else resolve(result!)
            }
        })
    })
}

function gzipAsync (
    f:Fflate,
    data:Uint8Array,
    opts:Record<string, unknown>
):Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
        f.gzip(data, opts, (err, result) => {
            if (err) reject(err)
            else {
                if (err !== null) reject(new Error('err not null'))
                else resolve(result!)
            }
        })
    })
}

function gunzipAsync (
    f:Fflate,
    data:Uint8Array
):Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
        f.gunzip(data, (err, result) => {
            if (err) reject(err)
            else {
                if (err !== null) reject(new Error('err not null'))
                else resolve(result!)
            }
        })
    })
}

function zlibAsync (
    f:Fflate,
    data:Uint8Array
):Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
        f.zlib(data, (err, result) => {
            if (err) reject(err)
            else {
                if (err !== null) reject(new Error('err not null'))
                else resolve(result!)
            }
        })
    })
}

function unzlibAsync (
    f:Fflate,
    data:Uint8Array
):Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
        f.unzlib(data, (err, result) => {
            if (err) reject(err)
            else {
                if (err !== null) reject(new Error('err not null'))
                else resolve(result!)
            }
        })
    })
}

function decompressAsync (
    f:Fflate,
    data:Uint8Array
):Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
        f.decompress(data, (err, result) => {
            if (err) reject(err)
            else {
                if (err !== null) reject(new Error('err not null'))
                else resolve(result!)
            }
        })
    })
}

export function asyncSuite (f:Fflate, label:string) {
    const fixtureList = fixtures()
    const fixtureNames = ['compressible', 'random', 'text', 'empty']

    for (let i = 0; i < fixtureList.length; i++) {
        const fixture = fixtureList[i]
        const name = fixtureNames[i]

        // Test deflate/inflate round-trip
        test(
      `${label} > deflate/inflate ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          try {
              const compressed = await withTimeout(
                  () => deflateAsync(f, orig)
              )
              const decompressed = await withTimeout(
                  () => inflateAsync(f, compressed)
              )
              t.ok(eq(decompressed, orig),
                  'deflate/inflate round-trip matches')
          } catch (e) {
              t.ok(false, `deflate/inflate failed: ${e}`)
          }
      }
        )

        // Test deflate sync/async match
        test(
      `${label} > deflate sync/async match ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          try {
              const syncDeflate = f.deflateSync(orig)
              const asyncDeflate = await withTimeout(
                  () => deflateAsync(f, orig)
              )
              t.ok(eq(asyncDeflate, syncDeflate),
                  'async deflate equals sync deflate')
          } catch (e) {
              t.ok(false, `deflate sync/async failed: ${e}`)
          }
      }
        )

        // Test gzip/gunzip round-trip
        test(
      `${label} > gzip/gunzip ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          try {
              const compressed = await withTimeout(
                  () => gzipAsync(f, orig, { mtime:MTIME })
              )
              const decompressed = await withTimeout(
                  () => gunzipAsync(f, compressed)
              )
              t.ok(eq(decompressed, orig),
                  'gzip/gunzip round-trip matches')
          } catch (e) {
              t.ok(false, `gzip/gunzip failed: ${e}`)
          }
      }
        )

        // Test gzip sync/async match
        test(
      `${label} > gzip sync/async match ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          try {
          // Pin the mtime: gzh writes Date.now()/1000 into header
          // bytes 4-7 (src/index.ts:1199), so an unpinned comparison
          // fails whenever the two calls straddle a second boundary.
          // Measured on the empty fixture: 4 mismatches in 200 at
          // offset 4 unpinned, 0 in 200 pinned.
              const syncGzip = f.gzipSync(orig, { mtime:MTIME })
              const asyncGzip = await withTimeout(
                  () => gzipAsync(f, orig, { mtime:MTIME })
              )
              const diff = firstDiff(asyncGzip, syncGzip)
              t.ok(diff === null, diff ?
            `async gzip differs: ${diff}` :
                  'async gzip equals sync gzip')
          } catch (e) {
              t.ok(false, `gzip sync/async failed: ${e}`)
          }
      }
        )

        // Test zlib/unzlib round-trip
        test(
      `${label} > zlib/unzlib ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          try {
              const compressed = await withTimeout(
                  () => zlibAsync(f, orig)
              )
              const decompressed = await withTimeout(
                  () => unzlibAsync(f, compressed)
              )
              t.ok(eq(decompressed, orig),
                  'zlib/unzlib round-trip matches')
          } catch (e) {
              t.ok(false, `zlib/unzlib failed: ${e}`)
          }
      }
        )

        // Test zlib sync/async match
        test(
      `${label} > zlib sync/async match ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          try {
              const syncZlib = f.zlibSync(orig)
              const asyncZlib = await withTimeout(
                  () => zlibAsync(f, orig)
              )
              t.ok(eq(asyncZlib, syncZlib),
                  'async zlib equals sync zlib')
          } catch (e) {
              t.ok(false, `zlib sync/async failed: ${e}`)
          }
      }
        )

        // Test decompress auto-detection
        if (name !== 'empty') {
            test(
        `${label} > decompress deflate ${name}`,
        async t => {
            const orig = new Uint8Array(fixture)
            try {
                const compressed = f.deflateSync(orig)
                const decompressed = await withTimeout(
                    () => decompressAsync(f, compressed)
                )
                t.ok(eq(decompressed, orig),
                    'decompress deflate works')
            } catch (e) {
                t.ok(false, `decompress deflate failed: ${e}`)
            }
        }
            )

            test(
        `${label} > decompress gzip ${name}`,
        async t => {
            const orig = new Uint8Array(fixture)
            try {
                const compressed = f.gzipSync(orig)
                const decompressed = await withTimeout(
                    () => decompressAsync(f, compressed)
                )
                t.ok(eq(decompressed, orig),
                    'decompress gzip works')
            } catch (e) {
                t.ok(false, `decompress gzip failed: ${e}`)
            }
        }
            )

            test(
        `${label} > decompress zlib ${name}`,
        async t => {
            const orig = new Uint8Array(fixture)
            try {
                const compressed = f.zlibSync(orig)
                const decompressed = await withTimeout(
                    () => decompressAsync(f, compressed)
                )
                t.ok(eq(decompressed, orig),
                    'decompress zlib works')
            } catch (e) {
                t.ok(false, `decompress zlib failed: ${e}`)
            }
        }
            )
        }
    }
}
