import { existsSync, readFile, writeFile } from 'fs'
import { resolve as resolvePath } from 'path'
import { fileURLToPath, pathToFileURL } from 'url'
import { get } from 'https'
import { test } from '@substrate-system/tapzero'
import type { Test } from '@substrate-system/tapzero'
import { performance } from 'perf_hooks'
import { Worker } from 'worker_threads'

// ESM has no __dirname
const here = fileURLToPath(new URL('.', import.meta.url))

const testFiles = {
    basic:Buffer.from('Hello world!'),
    text:'https://www.gutenberg.org/files/2701/old/moby10b.txt',
    smallImage:'https://hlevkin.com/hlevkin/TestImages/new/Rainier.bmp',
    image:'https://www.hlevkin.com/hlevkin/TestImages/new/Maltese.bmp',
    largeImage:'https://www.hlevkin.com/hlevkin/TestImages/new/Sunrise.bmp'
}

const testZipFiles = {
    model3D:'https://raw.githubusercontent.com/mrdoob/three.js/dev/examples/models/kmz/Box.kmz',
    largeModel3D:'https://raw.githubusercontent.com/mrdoob/three.js/dev/examples/models/3mf/truck.3mf'
}

const dlCached = async <T extends Record<string, string|Buffer>>(files:T) => {
    const res = {} as Record<keyof T, Buffer>
    for (const name in files) {
        let data:string|Buffer = files[name]
        if (typeof data === 'string') {
            const fn = resolvePath(here, 'data', name)
            if (!existsSync(fn)) {
                console.log('\nDownloading ' + data + '...')
                data = await new Promise((resolve, reject) => get(data as string, res => {
                    const len = +res.headers['content-length']
                    const buf = Buffer.allocUnsafe(len)
                    let i = 0
                    res.on('data', chunk => {
                        buf.set(chunk, i)
                        console.log((100 * (i += chunk.length) / len).toFixed(1) + '%\x1B[1A')
                    })
                    res.on('error', reject)
                    res.on('end', () => {
                        console.log('Complete')
                        writeFile(fn, buf, () => resolve(buf))
                    })
                }))
            } else {
                data = await new Promise((resolve, reject) =>
                    readFile(fn, (err, buf) => err ? reject(err) : resolve(buf))
                )
            }
        }
        res[name as keyof T] = data as Buffer
    }
    return res
}

const testFilesPromise = dlCached(testFiles)
const testZipFilesPromise = dlCached(testZipFiles)

export type TestHandler = (
    file:Buffer,
    name:string,
    resetTimer:() => void,
    t:Test
) => unknown|Promise<unknown>

export const testSuites = async <T extends Record<string, TestHandler>, D extends 'zip'|'default' = 'default'>(suites:T, type?:D) => {
    type DK = keyof (D extends 'zip' ? typeof testZipFiles : typeof testFiles)
    const tf = type === 'zip' ? testZipFiles : testFiles
    const tfp = type === 'zip' ? testZipFilesPromise : testFilesPromise
    const perf = {} as Record<keyof T, Promise<Record<DK, number>>>
    for (const k in suites) {
    // The executor must stay synchronous so that every test is registered
    // during module evaluation, before tapzero starts the run.
        perf[k] = new Promise(resolve => {
            const localPerf = {} as Record<DK, number>
            for (const name in tf) {
                test(k + ' > ' + name, async t => {
                    const localTestFiles = (await tfp) as unknown as Record<DK, Buffer>
                    let ts = performance.now()
                    await suites[k](localTestFiles[name as keyof typeof localTestFiles], name, () => {
                        ts = performance.now()
                    }, t)
                    localPerf[name as keyof typeof localPerf] = performance.now() - ts
                })
            }
            // Tests run in registration order, so this settles only once every
            // test in the suite above has finished.
            test(k + ' > collect timings', () => {
                resolve(localPerf)
            })
        })
    }
    const resolvedPerf = {} as Record<keyof T, Record<DK, number>>
    for (const k in suites) resolvedPerf[k] = await perf[k]
    return resolvedPerf
}

export const stream = (src:Uint8Array, dst:{
    push(dat:Uint8Array, final:boolean):void;
}) => {
    for (let i = 0; i < src.length;) {
        const off = Math.floor(Math.random() * Math.min(131072, src.length >>> 3))
        dst.push(src.slice(i, i + off), (i += off) >= src.length)
    }
}

// create worker string
//
// Workers stay CommonJS: they are spawned with { eval: true }, which
// Node evaluates as CJS. That is deliberate. A data: URL module worker
// cannot resolve bare specifiers (ERR_UNSUPPORTED_RESOLVE_REQUEST), and
// every comparison library below is loaded by bare specifier.
//
// fflate is now an ES module, so it is loaded with dynamic import() of
// an absolute file: URL. A CommonJS worker can do that.
const cws = (pkg:string, method:string = '_default') => {
    const load = pkg.startsWith('file:') ?
      `await import(${JSON.stringify(pkg)})` :
      `require(${JSON.stringify(pkg)})`

    // Normalises three shapes: an ESM namespace with named exports, a
    // CJS module whose export is the function, and an ESM default.
    const target = method === '_default' ?
        '(_m.default ?? _m)' :
      `(_m.default ?? _m).${method}`

    return `
    const { parentPort, workerData } = require('worker_threads');
    (async () => {
      try {
        const _m = ${load};
        const args = Array.isArray(workerData) ?
            workerData :
            [workerData];
        const buf = ${target}(...args);
        parentPort.postMessage(buf, [buf.buffer]);
      } catch (err) {
        const errPayload = err instanceof Error ?
            { name: err.name, message: err.message, stack: err.stack } :
            { name: 'Error', message: String(err) };
        parentPort.postMessage({ err: errPayload });
      }
    })();
  `
}

export type Workerized = (workerData:Uint8Array|[Uint8Array, object], transferable?:ArrayBuffer[]) => WorkerizedResult
export interface WorkerizedResult extends PromiseLike<Uint8Array<ArrayBuffer>> {
    timeout(ms:number):void;
};

// Worker creator
const wc = (pkg:string, method?:string):Workerized => {
    const str = cws(pkg, method)
    return (workerData, transferable) => {
        const worker = new Worker(str, {
            eval:true,
            workerData,
            transferList:transferable
        })
        let terminated = false
        return {
            timeout (ms:number) {
                const tm = setTimeout(() => {
                    worker.terminate()
                    terminated = true
                }, ms)
                worker.once('message', () => clearTimeout(tm))
            },
            then (res, rej) {
                return new Promise((resolve, reject) => {
                    worker
                        .once('message', msg => {
                            if (msg.err) {
                                return reject(Object.assign(
                                    new Error(msg.err.message),
                                    msg.err
                                ))
                            }
                            resolve(msg)
                        })
                        .once('error', reject)
                        .once('exit', code => {
                            if (terminated) reject(new Error('Timed out'))
                            else if (code !== 0) reject(new Error('Exited with status code ' + code))
                        })
                }).then(res, rej)
            }
        }
    }
}

// Workers load this by absolute URL, which bypasses the "exports" map.
// It is an ES module now, so cws() reaches it with dynamic import().
const fflate = pathToFileURL(
    resolvePath(here, '..', 'dist', 'node', 'index.js')
).href

export const workers = {
    fflate:{
        deflate:wc(fflate, 'deflateSync'),
        inflate:wc(fflate, 'inflateSync'),
        gzip:wc(fflate, 'gzipSync'),
        gunzip:wc(fflate, 'gunzipSync'),
        zlib:wc(fflate, 'zlibSync'),
        unzlib:wc(fflate, 'unzlibSync'),
        // zip and unzip are not exercised by the node harness. The skip at
        // test/2-perf.ts:26 is pre-existing upstream and unchanged by this
        // branch.
        //
        // unzip additionally cannot round trip through wc(): it returns an
        // Unzipped plain object with no .buffer, so this harness's fixed
        // postMessage(buf, [buf.buffer]) transfer list rejects it with
        // DataCloneError.
        //
        // ZIP coverage lands in Phase 6's browser suite.
        zip:wc(fflate, 'zipSync'),
        unzip:wc(fflate, 'unzipSync')
    },
    pako:{
        deflate:wc('pako', 'deflateRaw'),
        inflate:wc('pako', 'inflateRaw'),
        gzip:wc('pako', 'gzip'),
        gunzip:wc('pako', 'ungzip'),
        zlib:wc('pako', 'deflate'),
        unzlib:wc('pako', 'inflate')
    },
    uzip:{
        deflate:wc('uzip', 'deflateRaw'),
        inflate:wc('uzip', 'inflateRaw')
    },
    tinyInflate:{
        inflate:wc('tiny-inflate')
    },
    zlib:{
        deflate:wc('zlib', 'deflateRawSync'),
        inflate:wc('zlib', 'inflateRawSync'),
        gzip:wc('zlib', 'gzipSync'),
        gunzip:wc('zlib', 'gunzipSync'),
        zlib:wc('zlib', 'deflateSync'),
        unzlib:wc('zlib', 'inflateSync')
    }
}

export const bClone = (buf:Buffer) => {
    const clone = Buffer.allocUnsafeSlow(buf.length)
    clone.set(buf)
    return clone
}
