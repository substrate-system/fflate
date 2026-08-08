import type * as fflate from '../../src/index.js'

export type Fflate = typeof fflate

// Deterministic pseudo-random source shared by fixtures() and
// largeIncompressible(). Returns a float in [0, 1).
//
// The generator state is xorshift128, kept entirely in bitwise 32-bit
// operations. That matters: the obvious LCG,
// `state = (state * 1103515245 + 12345) & 0x7fffffff`, silently
// degenerates in JS because the product exceeds 2**53 and the low bits
// are lost to float rounding before the mask applies. An earlier
// version of largeIncompressible() used it and deflated 600000 bytes
// to 19536, which sent every caller down the synchronous branch. Do
// not reintroduce an LCG here.
//
// Only the division at the end leaves integer arithmetic, and it is
// exact: 0x100000000 is a power of two, so callers doing
// `Math.floor(gen() * 256)` get the top 8 bits of the state with no
// rounding.
function xorshift128 ():() => number {
    let x = 123456789
    let y = 362436069
    let z = 521288629
    let w = 88675123
    return () => {
        const t = x ^ (x << 11)
        x = y; y = z; z = w
        w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8))
        return (w >>> 0) / 0x100000000
    }
}

export function fixtures ():Uint8Array[] {
    // Highly compressible: long run of one byte
    const compressible = new Uint8Array(1000)
    compressible.fill(65) // 'A'

    // Incompressible: deterministic pseudo-random sequence. Measured:
    // 512 bytes -> 517 deflated, 222 distinct byte values.
    const random = new Uint8Array(512)
    const xs = xorshift128()
    for (let i = 0; i < random.length; i++) {
        random[i] = Math.floor(xs() * 256)
    }

    // Text: UTF-8 including multi-byte characters
    const text = 'Hello, World! 你好世界 🌍 Привет мир'
    const textBytes = new Uint8Array(
        new TextEncoder().encode(text)
    )

    // Empty
    const empty = new Uint8Array(0)

    return [compressible, random, textBytes, empty]
}

// Large fixtures for async worker path tests.
//
// The two async paths want fixtures with OPPOSITE properties, so they
// get separate helpers. Measured sizes are stated below; if either
// helper changes, re-measure rather than assuming.
//
//   zip()    offloads to a worker when originalSize >= 160000
//   unzip()  offloads to a worker when originalSize >= 524288 AND
//            compressedSize <= 0.8 * originalSize  -> COMPRESSIBLE
//   AsyncUnzipInflate picks AsyncInflate when the COMPRESSED size
//            passed as `sz` is >= 320000            -> INCOMPRESSIBLE

// 600000 bytes that deflate to 599. Clears both the 160000 zip
// threshold and the 524288 unzip threshold, and 599 is far under
// 0.8 * 600000, so unzip() takes its worker branch.
export function largeCompressible ():Uint8Array {
    const fixture = new Uint8Array(600000)
    fixture.fill(66) // 'B'
    return fixture
}

// 600000 bytes that deflate to 600125, i.e. genuinely incompressible,
// so AsyncUnzipInflate receives sz >= 320000 and takes its async
// branch. See xorshift128() above for why the generator matters here,
// and re-measure with deflateSync if you change this.
export function largeIncompressible ():Uint8Array {
    const xs = xorshift128()
    const fixture = new Uint8Array(600000)
    for (let i = 0; i < fixture.length; i++) {
        fixture[i] = Math.floor(xs() * 256)
    }
    return fixture
}

// Joins the chunks a streaming class emits through ondata. Streaming
// compressors emit many chunks; keeping only the last one and comparing
// it against one-shot output is a silent false failure.
export function concat (parts:Uint8Array[]):Uint8Array {
    let total = 0
    for (const p of parts) total += p.length
    const out = new Uint8Array(total)
    let offset = 0
    for (const p of parts) {
        out.set(p, offset)
        offset += p.length
    }
    return out
}

export function eq (a:Uint8Array, b:Uint8Array):boolean {
    if (a.length !== b.length) return false
    for (let i = 0; i < a.length; i++) {
        if (a[i] !== b[i]) return false
    }
    return true
}

export function chunks (
    data:Uint8Array,
    n:number
):Uint8Array[] {
    const result:Uint8Array[] = []
    const chunkSize = Math.ceil(data.length / n)
    for (let i = 0; i < data.length; i += chunkSize) {
        result.push(data.slice(i, i + chunkSize))
    }
    return result
}

// Reports the first byte offset where two Uint8Arrays differ, or null if
// equal. On mismatch, includes lengths and byte values for debugging.
export function firstDiff (
    a:Uint8Array,
    b:Uint8Array
):string|null {
    if (a.length !== b.length) {
        return `length mismatch: ${a.length} vs ${b.length}`
    }
    for (let i = 0; i < a.length; i++) {
        if (a[i] !== b[i]) {
            return (
        `byte mismatch at offset ${i}: ` +
        `got 0x${a[i].toString(16).padStart(2, '0')} ` +
        `expected 0x${b[i].toString(16).padStart(2, '0')}`
            )
        }
    }
    return null
}

// Race a function against a timeout, returning the result or rejecting
// if the timeout fires first.
//
// The default MUST stay below tapout's auto-finish window, or a stalled
// test truncates the run instead of reporting. tapout ends the run
// after `max(500, min(3000, floor(timeout * 0.2)))` ms with no console
// output (node_modules/@substrate-system/tapout/dist/test-harness.js:29)
// and resets that timer on every line printed. Our suites are silent while
// polling, so a hang goes quiet and auto-finish fires -- dropping every
// later test, including the entire minified pass, which is the one
// thing this suite exists to run. `test:browser` therefore passes
// `--timeout 30000`, putting the window at its 3000 ms cap, and this
// default sits under it so the rejection lands first: the catch prints
// `not ok`, which resets the window, and the run carries on.
//
// tapout's default timeout of 5000 would put the window at 1000 ms,
// below this value, and restore the truncating behaviour.
export const DEFAULT_TIMEOUT_MS = 2000

// Slack demanded on top of DEFAULT_TIMEOUT_MS by harness.ts. The two
// timers do not start together: auto-finish starts at the previous
// test's output line, withTimeout only once the test body reaches its
// await, after any setup work.
export const TIMEOUT_MARGIN_MS = 500

// Mirrors tapout's own formula at
// node_modules/@substrate-system/tapout/dist/test-harness.js:29.
// `timeoutMs` is whatever `--timeout` the test:browser script passes.
export function autoFinishWindow (timeoutMs:number):number {
    return Math.max(500, Math.min(3000, Math.floor(timeoutMs * 0.2)))
}

export function withTimeout<T> (
    fn:() => Promise<T>,
    ms:number = DEFAULT_TIMEOUT_MS
):Promise<T> {
    let timerId:ReturnType<typeof setTimeout>|null = null
    return Promise.race([
        fn(),
        new Promise<T>((_resolve, reject) => {
            timerId = setTimeout(
                () => reject(new Error('timeout')),
                ms
            )
        })
    ]).finally(() => {
        if (timerId !== null) {
            clearTimeout(timerId)
        }
    })
}
