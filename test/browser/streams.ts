import { test } from '@substrate-system/tapzero'
import type { Fflate } from './util.js'
import {
    fixtures, eq, chunks, concat, withTimeout
} from './util.js'

// gzip headers carry a 1-second-resolution mtime, so any byte-identity
// comparison must pin it on both sides. See the AsyncGzip test.
const MTIME = new Date('2020-01-01T00:00:00Z')

export function streamSuite (f:Fflate, label:string) {
    const fixtureList = fixtures()
    const fixtureNames = ['compressible', 'random', 'text', 'empty']

    for (let i = 0; i < fixtureList.length; i++) {
        const fixture = fixtureList[i]
        const name = fixtureNames[i]

        // Sync Deflate
        test(`${label} > Deflate sync ${name}`, async t => {
            const orig = new Uint8Array(fixture)
            const compressed:Uint8Array[] = []

            const deflater = new f.Deflate((chunk) => {
                compressed.push(new Uint8Array(chunk))
            })

            const numChunks = orig.length > 0 ? 3 : 1
            const chunkList = chunks(orig, numChunks)
            for (let j = 0; j < chunkList.length; j++) {
                const isLast = j === chunkList.length - 1
                deflater.push(chunkList[j], isLast)
            }
            if (orig.length === 0) {
                deflater.push(new Uint8Array(0), true)
            }

            const allCompressed = concat(compressed)

            const decompressed = f.inflateSync(allCompressed)
            t.ok(eq(decompressed, orig),
                'sync Deflate round-trip matches original')
        })

        // Async AsyncDeflate vs sync Deflate
        test(
      `${label} > AsyncDeflate vs sync ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)

          // Sync reference
          const syncCompressed:Uint8Array[] = []
          const syncDeflater = new f.Deflate((chunk) => {
              syncCompressed.push(new Uint8Array(chunk))
          })
          const numChunks = orig.length > 0 ? 3 : 1
          const chunkList = chunks(orig, numChunks)
          for (let j = 0; j < chunkList.length; j++) {
              const isLast = j === chunkList.length - 1
              syncDeflater.push(chunkList[j], isLast)
          }
          if (orig.length === 0) {
              syncDeflater.push(new Uint8Array(0), true)
          }

          const syncAll = concat(syncCompressed)

          // Async version
          try {
              const asyncCompressed:Uint8Array[] = []
              let deflateErr:Error|null = null
              let deflateSettled = false
              const asyncDeflater = new f.AsyncDeflate(
                  (err, chunk, final) => {
                      if (err) {
                          deflateErr = err
                          deflateSettled = true
                          return
                      }
                      asyncCompressed.push(new Uint8Array(chunk))
                      if (final) deflateSettled = true
                  }
              )

              for (let j = 0; j < chunkList.length; j++) {
                  const isLast = j === chunkList.length - 1
                  asyncDeflater.push(chunkList[j], isLast)
              }
              if (orig.length === 0) {
                  asyncDeflater.push(new Uint8Array(0), true)
              }

              await withTimeout(() =>
                  new Promise<void>(resolve => {
                      const check = () => {
                          if (deflateSettled) {
                              resolve()
                          } else {
                              setTimeout(check, 10)
                          }
                      }
                      check()
                  })
              )

              t.ok(deflateErr === null, deflateErr ?
            `AsyncDeflate callback error: ${deflateErr}` :
                  'AsyncDeflate reported no error')

              const asyncAll = concat(asyncCompressed)

              t.ok(eq(asyncAll, syncAll),
                  'async AsyncDeflate equals sync Deflate')
          } catch (e) {
              t.ok(false, `AsyncDeflate failed: ${e}`)
          }
      }
        )

        // Sync Inflate
        test(`${label} > Inflate sync ${name}`, async t => {
            const orig = new Uint8Array(fixture)
            const compressed = f.deflateSync(orig)

            const decompressed:Uint8Array[] = []
            const inflater = new f.Inflate((chunk) => {
                decompressed.push(new Uint8Array(chunk))
            })

            const compChunks = chunks(compressed, 2)
            for (let j = 0; j < compChunks.length; j++) {
                const isLast = j === compChunks.length - 1
                inflater.push(compChunks[j], isLast)
            }

            const result = concat(decompressed)

            t.ok(eq(result, orig),
                'sync Inflate round-trip matches original')
        })

        // Async AsyncInflate vs original
        test(
      `${label} > AsyncInflate equals original ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          const compressed = f.deflateSync(orig)

          try {
              const asyncDecompressed:Uint8Array[] = []
              let inflateErr:Error|null = null
              let inflateSettled = false
              const asyncInflater = new f.AsyncInflate(
                  (err, chunk, final) => {
                      if (err) {
                          inflateErr = err
                          inflateSettled = true
                          return
                      }
                      asyncDecompressed.push(new Uint8Array(chunk))
                      if (final) inflateSettled = true
                  }
              )

              const compChunks = chunks(compressed, 2)
              for (let j = 0; j < compChunks.length; j++) {
                  const isLast = j === compChunks.length - 1
                  asyncInflater.push(compChunks[j], isLast)
              }

              await withTimeout(() =>
                  new Promise<void>(resolve => {
                      const check = () => {
                          if (inflateSettled) {
                              resolve()
                          } else {
                              setTimeout(check, 10)
                          }
                      }
                      check()
                  })
              )

              t.ok(inflateErr === null, inflateErr ?
            `AsyncInflate callback error: ${inflateErr}` :
                  'AsyncInflate reported no error')

              const asyncResult = concat(asyncDecompressed)

              t.ok(eq(asyncResult, orig),
                  'async AsyncInflate equals sync Inflate')
          } catch (e) {
              t.ok(false, `AsyncInflate failed: ${e}`)
          }
      }
        )

        // Sync Gzip
        test(`${label} > Gzip sync ${name}`, async t => {
            const orig = new Uint8Array(fixture)
            const gzipped:Uint8Array[] = []

            const gzipper = new f.Gzip((chunk) => {
                gzipped.push(new Uint8Array(chunk))
            })

            const numChunks = orig.length > 0 ? 3 : 1
            const chunkList = chunks(orig, numChunks)
            for (let j = 0; j < chunkList.length; j++) {
                const isLast = j === chunkList.length - 1
                gzipper.push(chunkList[j], isLast)
            }
            if (orig.length === 0) {
                gzipper.push(new Uint8Array(0), true)
            }

            const allGzipped = concat(gzipped)

            const decompressed = f.gunzipSync(allGzipped)
            t.ok(eq(decompressed, orig),
                'sync Gzip round-trip matches original')
        })

        // Async AsyncGzip vs sync Gzip
        test(
      `${label} > AsyncGzip vs sync ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          const numChunks = orig.length > 0 ? 3 : 1
          const chunkList = chunks(orig, numChunks)

          // Compare against the SYNC STREAM over the same chunking,
          // not against one-shot gzipSync. Chunked streaming emits
          // different block boundaries than a single-shot compression:
          // on incompressible input the sync stream itself differs from
          // the one-shot output (561 vs 535 bytes, measured), so
          // one-shot equality was never a real invariant. Async equals
          // sync at identical chunking is, and it is the comparison
          // that isolates the worker path.
          try {
          // Under withTimeout and inside the try: a throw from
          // s.push, or a stream that never emits final, would
          // otherwise reject or hang outside any handler, and
          // tapzero aborts the whole run on an escaped rejection --
          // taking the minified pass with it.
          //
          // Pin the mtime on both sides. gzh writes Date.now()/1000
          // into header bytes 4-7 (src/index.ts:1199), so an unpinned
          // byte comparison fails whenever the two streams straddle a
          // second boundary. Measured one-shot: 4 mismatches in 200
          // at offset 4 unpinned, 0 in 200 pinned.
              const syncGzip = await withTimeout(() =>
                  new Promise<Uint8Array>(resolve => {
                      const syncParts:Uint8Array[] = []
                      const s = new f.Gzip({ mtime:MTIME }, (chunk, final) => {
                          syncParts.push(new Uint8Array(chunk))
                          if (final) resolve(concat(syncParts))
                      })
                      if (orig.length === 0) s.push(new Uint8Array(0), true)
                      for (let j = 0; j < chunkList.length; j++) {
                          s.push(chunkList[j], j === chunkList.length - 1)
                      }
                  })
              )

              const parts:Uint8Array[] = []
              let asyncGzipped:Uint8Array|null = null
              let gzipError:Error|null = null
              let gzipSettled = false
              const gzipper = new f.AsyncGzip(
                  { mtime:MTIME },
                  (err, chunk, final) => {
                      if (err) {
                          gzipError = err
                          gzipSettled = true
                          return
                      }
                      parts.push(new Uint8Array(chunk))
                      // Accumulate: a streaming compressor emits many
                      // chunks, and only the full join equals one-shot
                      // output. Settle only once final arrives.
                      if (final) {
                          asyncGzipped = concat(parts)
                          gzipSettled = true
                      }
                  }
              )

              for (let j = 0; j < chunkList.length; j++) {
                  const isLast = j === chunkList.length - 1
                  gzipper.push(chunkList[j], isLast)
              }
              if (orig.length === 0) {
                  gzipper.push(new Uint8Array(0), true)
              }

              // Poll on gzipSettled, which the error branch sets too.
              // Polling on asyncGzipped instead would make both
              // assertions below unfailable: the poll would only resolve
              // once the thing being asserted was already true, and an
              // error would spin until the timeout.
              await withTimeout(() =>
                  new Promise<void>((resolve) => {
                      const check = () => {
                          if (gzipSettled) {
                              resolve()
                          } else {
                              setTimeout(check, 10)
                          }
                      }
                      check()
                  })
              )

              t.ok(gzipError === null, gzipError ?
            `AsyncGzip callback error: ${gzipError}` :
                  'AsyncGzip reported no error')
              t.ok(asyncGzipped !== null && eq(asyncGzipped, syncGzip),
                  'async AsyncGzip equals sync Gzip')
          } catch (e) {
              t.ok(false, `AsyncGzip failed: ${e}`)
          }
      }
        )

        // Sync Gunzip
        test(`${label} > Gunzip sync ${name}`, async t => {
            const orig = new Uint8Array(fixture)
            const gzipped = f.gzipSync(orig)

            const gunzipped:Uint8Array[] = []
            const gunzipper = new f.Gunzip((chunk) => {
                gunzipped.push(new Uint8Array(chunk))
            })

            const gzipChunks = chunks(gzipped, 2)
            for (let j = 0; j < gzipChunks.length; j++) {
                const isLast = j === gzipChunks.length - 1
                gunzipper.push(gzipChunks[j], isLast)
            }

            const result = concat(gunzipped)

            t.ok(eq(result, orig),
                'sync Gunzip round-trip matches original')
        })

        // Async AsyncGunzip vs original
        test(
      `${label} > AsyncGunzip equals original ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          const gzipped = f.gzipSync(orig)

          try {
              const asyncGunzipped:Uint8Array[] = []
              let gunzipErr:Error|null = null
              let gunzipSettled = false
              const gunzipper = new f.AsyncGunzip(
                  (err, chunk, final) => {
                      if (err) {
                          gunzipErr = err
                          gunzipSettled = true
                          return
                      }
                      asyncGunzipped.push(new Uint8Array(chunk))
                      if (final) gunzipSettled = true
                  }
              )

              const gzipChunks = chunks(gzipped, 2)
              for (let j = 0; j < gzipChunks.length; j++) {
                  const isLast = j === gzipChunks.length - 1
                  gunzipper.push(gzipChunks[j], isLast)
              }

              await withTimeout(() =>
                  new Promise<void>(resolve => {
                      const check = () => {
                          if (gunzipSettled) {
                              resolve()
                          } else {
                              setTimeout(check, 10)
                          }
                      }
                      check()
                  })
              )

              t.ok(gunzipErr === null, gunzipErr ?
            `AsyncGunzip callback error: ${gunzipErr}` :
                  'AsyncGunzip reported no error')

              const asyncResult = concat(asyncGunzipped)

              t.ok(eq(asyncResult, orig),
                  'async AsyncGunzip equals original')
          } catch (e) {
              t.ok(false, `AsyncGunzip failed: ${e}`)
          }
      }
        )

        // Sync Zlib
        test(`${label} > Zlib sync ${name}`, async t => {
            const orig = new Uint8Array(fixture)
            const zlibbed:Uint8Array[] = []

            const zliber = new f.Zlib((chunk) => {
                zlibbed.push(new Uint8Array(chunk))
            })

            const numChunks = orig.length > 0 ? 3 : 1
            const chunkList = chunks(orig, numChunks)
            for (let j = 0; j < chunkList.length; j++) {
                const isLast = j === chunkList.length - 1
                zliber.push(chunkList[j], isLast)
            }
            if (orig.length === 0) {
                zliber.push(new Uint8Array(0), true)
            }

            const allZlibbed = concat(zlibbed)

            const decompressed = f.unzlibSync(allZlibbed)
            t.ok(eq(decompressed, orig),
                'sync Zlib round-trip matches original')
        })

        // Async AsyncZlib vs sync Zlib
        test(
      `${label} > AsyncZlib vs sync ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          const numChunks = orig.length > 0 ? 3 : 1
          const chunkList = chunks(orig, numChunks)

          // Compare against the SYNC STREAM over the same chunking,
          // not against one-shot zlibSync. Chunked streaming emits
          // different block boundaries than a single-shot compression:
          // on incompressible input the sync stream itself differs from
          // the one-shot output (561 vs 535 bytes, measured), so
          // one-shot equality was never a real invariant. Async equals
          // sync at identical chunking is, and it is the comparison
          // that isolates the worker path.
          try {
          // Under withTimeout and inside the try: a throw from
          // s.push, or a stream that never emits final, would
          // otherwise reject or hang outside any handler, and
          // tapzero aborts the whole run on an escaped rejection --
          // taking the minified pass with it.
              const syncZlib = await withTimeout(() =>
                  new Promise<Uint8Array>(resolve => {
                      const syncParts:Uint8Array[] = []
                      const s = new f.Zlib((chunk, final) => {
                          syncParts.push(new Uint8Array(chunk))
                          if (final) resolve(concat(syncParts))
                      })
                      if (orig.length === 0) s.push(new Uint8Array(0), true)
                      for (let j = 0; j < chunkList.length; j++) {
                          s.push(chunkList[j], j === chunkList.length - 1)
                      }
                  })
              )

              const parts:Uint8Array[] = []
              let asyncZlibbed:Uint8Array|null = null
              let zlibError:Error|null = null
              let zlibSettled = false
              const zliber = new f.AsyncZlib((err, chunk, final) => {
                  if (err) {
                      zlibError = err
                      zlibSettled = true
                      return
                  }
                  parts.push(new Uint8Array(chunk))
                  // Accumulate: a streaming compressor emits many
                  // chunks, and only the full join equals one-shot
                  // output. Settle only once final arrives.
                  if (final) {
                      asyncZlibbed = concat(parts)
                      zlibSettled = true
                  }
              })

              for (let j = 0; j < chunkList.length; j++) {
                  const isLast = j === chunkList.length - 1
                  zliber.push(chunkList[j], isLast)
              }
              if (orig.length === 0) {
                  zliber.push(new Uint8Array(0), true)
              }

              // Poll on zlibSettled, not on asyncZlibbed -- see the
              // AsyncGzip test above for why.
              await withTimeout(() =>
                  new Promise<void>((resolve) => {
                      const check = () => {
                          if (zlibSettled) {
                              resolve()
                          } else {
                              setTimeout(check, 10)
                          }
                      }
                      check()
                  })
              )

              t.ok(zlibError === null, zlibError ?
            `AsyncZlib callback error: ${zlibError}` :
                  'AsyncZlib reported no error')
              t.ok(asyncZlibbed !== null && eq(asyncZlibbed, syncZlib),
                  'async AsyncZlib equals sync Zlib')
          } catch (e) {
              t.ok(false, `AsyncZlib failed: ${e}`)
          }
      }
        )

        // Sync Unzlib
        test(`${label} > Unzlib sync ${name}`, async t => {
            const orig = new Uint8Array(fixture)
            const zlibbed = f.zlibSync(orig)

            const unzlibbed:Uint8Array[] = []
            const unzliber = new f.Unzlib((chunk) => {
                unzlibbed.push(new Uint8Array(chunk))
            })

            const zlibChunks = chunks(zlibbed, 2)
            for (let j = 0; j < zlibChunks.length; j++) {
                const isLast = j === zlibChunks.length - 1
                unzliber.push(zlibChunks[j], isLast)
            }

            const result = concat(unzlibbed)

            t.ok(eq(result, orig),
                'sync Unzlib round-trip matches original')
        })

        // Async AsyncUnzlib vs original
        test(
      `${label} > AsyncUnzlib equals original ${name}`,
      async t => {
          const orig = new Uint8Array(fixture)
          const zlibbed = f.zlibSync(orig)

          try {
              const asyncUnzlibbed:Uint8Array[] = []
              let unzlibErr:Error|null = null
              let unzlibSettled = false
              const unzliber = new f.AsyncUnzlib(
                  (err, chunk, final) => {
                      if (err) {
                          unzlibErr = err
                          unzlibSettled = true
                          return
                      }
                      asyncUnzlibbed.push(new Uint8Array(chunk))
                      if (final) unzlibSettled = true
                  }
              )

              const zlibChunks = chunks(zlibbed, 2)
              for (let j = 0; j < zlibChunks.length; j++) {
                  const isLast = j === zlibChunks.length - 1
                  unzliber.push(zlibChunks[j], isLast)
              }

              await withTimeout(() =>
                  new Promise<void>(resolve => {
                      const check = () => {
                          if (unzlibSettled) {
                              resolve()
                          } else {
                              setTimeout(check, 10)
                          }
                      }
                      check()
                  })
              )

              t.ok(unzlibErr === null, unzlibErr ?
            `AsyncUnzlib callback error: ${unzlibErr}` :
                  'AsyncUnzlib reported no error')

              const asyncResult = concat(asyncUnzlibbed)

              t.ok(eq(asyncResult, orig),
                  'async AsyncUnzlib equals original')
          } catch (e) {
              t.ok(false, `AsyncUnzlib failed: ${e}`)
          }
      }
        )

        // Sync Decompress
        test(`${label} > Decompress sync ${name}`, async t => {
            const orig = new Uint8Array(fixture)
            const deflated = f.deflateSync(orig)

            const decompressed:Uint8Array[] = []
            const decomposer = new f.Decompress((chunk) => {
                decompressed.push(new Uint8Array(chunk))
            })

            const defChunks = chunks(deflated, 2)
            for (let j = 0; j < defChunks.length; j++) {
                const isLast = j === defChunks.length - 1
                decomposer.push(defChunks[j], isLast)
            }

            const result = concat(decompressed)

            t.ok(eq(result, orig),
                'sync Decompress round-trip matches original')
        })

        // Empty final chunk test for async streaming
        if (name === 'empty') {
            test(
        `${label} > AsyncDeflate empty final chunk`,
        async t => {
            try {
                const dataCollected:Uint8Array[] = []
                let emptyFinalErr:Error|null = null
                let emptyFinalSettled = false
                const deflater = new f.AsyncDeflate(
                    (err, chunk, final) => {
                        if (err) {
                            emptyFinalErr = err
                            emptyFinalSettled = true
                            return
                        }
                        dataCollected.push(new Uint8Array(chunk))
                        if (final) emptyFinalSettled = true
                    }
                )

                // Push empty chunk with final flag
                deflater.push(new Uint8Array(0), true)

                await withTimeout(() =>
                    new Promise<void>(resolve => {
                        const check = () => {
                            if (emptyFinalSettled) {
                                resolve()
                            } else {
                                setTimeout(check, 10)
                            }
                        }
                        check()
                    })
                )

                t.ok(emptyFinalErr === null, emptyFinalErr ?
              `AsyncDeflate empty final callback error: ${emptyFinalErr}` :
                    'AsyncDeflate empty final reported no error')

                // Assert that the output is a valid DEFLATE stream that
                // inflates back to empty
                const compressed = concat(dataCollected)
                const inflated = f.inflateSync(compressed)
                t.ok(inflated.length === 0,
                    'AsyncDeflate empty final chunk inflates to empty')
            } catch (e) {
                t.ok(false, `empty final chunk failed: ${e}`)
            }
        }
            )
        }
    }

    // Compressed bytes split across chunk boundary (tests AsyncInflate
    // chunking, not character encoding state)
    test(
    `${label} > AsyncInflate compressed split across chunks`,
    async t => {
        // The emoji 🌍 is 4 bytes in UTF-8: F0 9F 8C 8D
        const emojiBytes = new Uint8Array([0xF0, 0x9F, 0x8C, 0x8D])

        const decompressed:Uint8Array[] = []
        let emojiErr:Error|null = null
        let emojiSettled = false
        const inflater = new f.AsyncInflate(
            (err, chunk, final) => {
                if (err) {
                    emojiErr = err
                    emojiSettled = true
                    return
                }
                decompressed.push(new Uint8Array(chunk))
                if (final) emojiSettled = true
            }
        )

        // First, compress the full emoji
        const compressed = f.deflateSync(emojiBytes)

        // Now decompress it in chunks that will split COMPRESSED bytes,
        // testing AsyncInflate's ability to handle fragmented input
        const chunkSize = Math.ceil(compressed.length / 2)
        for (let i = 0; i < compressed.length; i += chunkSize) {
            const end = Math.min(i + chunkSize, compressed.length)
            const isLast = end === compressed.length
            inflater.push(compressed.slice(i, end), isLast)
        }

        try {
            await withTimeout(() =>
                new Promise<void>(resolve => {
                    const check = () => {
                        if (emojiSettled) {
                            resolve()
                        } else {
                            setTimeout(check, 10)
                        }
                    }
                    check()
                })
            )

            t.ok(emojiErr === null, emojiErr ?
          `AsyncInflate callback error: ${emojiErr}` :
                'AsyncInflate reported no error')

            const result = concat(decompressed)

            t.ok(eq(result, emojiBytes),
                'emoji split across compressed chunks preserved')
        } catch (e) {
            t.ok(false, `AsyncInflate compressed split across chunks failed: ${e}`)
        }
    }
    )

    // EncodeUTF8 and DecodeUTF8: multi-byte character split across
    // DECODED stream chunk boundary. These classes hold state across
    // chunks to handle partial characters.
    test(
    `${label} > EncodeUTF8 multi-byte character`,
    async t => {
        // Multi-byte UTF-8 characters that will be split across chunks
        const text = 'Hello🌍'
        const encoded = f.strToU8(text)

        // Use EncodeUTF8 and push chunks that split the emoji
        const encoder = new f.EncodeUTF8()
        const chunksOut:Uint8Array[] = []

        encoder.ondata = (chunk) => {
            chunksOut.push(new Uint8Array(chunk))
        }

        // Push text in chunks: 'Hello' then '🌍'
        // 'Hello' is 5 bytes, emoji is 4 bytes
        encoder.push('Hello', false)
        encoder.push('🌍', true)

        const result = concat(chunksOut)

        t.ok(eq(result, encoded),
            'EncodeUTF8 multi-byte characters work')
    }
    )

    // DecodeUTF8: decode with multi-byte character split
    test(
    `${label} > DecodeUTF8 multi-byte character`,
    async t => {
        const text = 'Hello🌍'
        const encoded = f.strToU8(text)

        // Use DecodeUTF8 and push chunks that split the emoji
        // The emoji 🌍 is F0 9F 8C 8D
        // Split it: first chunk includes 'Hello' and F0 9F,
        // second chunk gets 8C 8D
        const decoder = new f.DecodeUTF8()
        const chunksOut:string[] = []

        decoder.ondata = (chunk) => {
            chunksOut.push(chunk)
        }

        // 'Hello' is 5 bytes, first 2 bytes of emoji
        const chunk1 = encoded.slice(0, 7)
        // Last 2 bytes of emoji
        const chunk2 = encoded.slice(7)

        decoder.push(chunk1, false)
        decoder.push(chunk2, true)

        const result = chunksOut.join('')
        t.ok(result === text,
            'DecodeUTF8 multi-byte character split works')
    }
    )
}
