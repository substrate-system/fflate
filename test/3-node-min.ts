import { test } from '@substrate-system/tapzero'

// Import both plain and minified node bundles
// We use dynamic imports and await to ensure they load
const importBundles = async () => {
    const plain = await import('../dist/node/index.js')
    // @ts-expect-error no declarations for minified bundle
    const min:typeof plain = await import('../dist/node/index.min.js')
    return { plain, min }
}

test('node bundle min: exports match plain', async t => {
    const { plain, min } = await importBundles()

    const plainKeys = Object.keys(plain).sort()
    const minKeys = Object.keys(min).sort()

    // The length guard is not redundant: without it two empty namespaces
    // would satisfy the comparison.
    const onlyPlain = plainKeys.filter(k => !minKeys.includes(k))
    const onlyMin = minKeys.filter(k => !plainKeys.includes(k))

    t.ok(
        plainKeys.length > 0 && !onlyPlain.length && !onlyMin.length,
    `export names match: plain ${plainKeys.length}, min ${minKeys.length}` +
      (onlyPlain.length ? `, only in plain: ${onlyPlain.join()}` : '') +
      (onlyMin.length ? `, only in min: ${onlyMin.join()}` : '')
    )
})

test('node bundle min: sync roundtrip', async t => {
    const { min } = await importBundles()

    const input = new TextEncoder().encode('hello sync world')
    const compressed = min.deflateSync(input)
    const decompressed = min.inflateSync(compressed)

    const match = input.length === decompressed.length &&
    input.every((v, i) => v === decompressed[i])

    t.ok(match, 'sync roundtrip preserves data')
})

test('node bundle min: async roundtrip large', async t => {
    const { min } = await importBundles()

    // Async deflate/inflate reach the worker unconditionally: deflate()
    // calls cbify() with no size gate (src/index.ts:1472). The size gates
    // at 160000 / 320000 / 524288 belong to the zip and unzip paths, not
    // this one. So the fixture's job is not to clear a threshold -- it is
    // to make the payload big and incompressible enough that the worker
    // does real DEFLATE work rather than a trivial run-length encode.
    //
    // Generate it with xorshift128 in exact 32-bit arithmetic. An LCG
    // degenerates here: the product exceeds 2**53 and the low bits are
    // lost to float rounding before the mask applies, yielding data that
    // deflates 30x and defeats the point of the fixture.
    const state = [123456789, 362436069, 521288629, 88675123]
    const fixture = new Uint8Array(600000)
    for (let i = 0; i < fixture.length; ++i) {
        const a = state[0]
        const d = state[3]
        const tn = a ^ (a << 11)
        state[0] = state[1]
        state[1] = state[2]
        state[2] = state[3]
        state[3] = d ^ (d >>> 19) ^ (tn ^ (tn >>> 8))
        fixture[i] = (state[3] >>> 0) & 0xff
    }

    // Assert the fixture's own property, so a drifting generator fails
    // loudly here instead of silently weakening the test below.
    const ratio = min.deflateSync(fixture).length / fixture.length
    t.ok(ratio > 0.95, `fixture is incompressible: ratio ${ratio}`)

    // Now test async. Hold on to the AsyncTerminable each call returns:
    // cbify only calls w.terminate() from inside the worker's own
    // callback (src/index.ts:1123), so a worker that never calls back --
    // the exact wcln failure this test guards -- keeps the node event
    // loop alive forever. Without terminating it here the suite reports
    // its failures, prints its plan line, and then hangs instead of
    // exiting, which CI sees as a job timeout rather than exit 1.
    let result:Uint8Array|null = null
    let error:Error|null = null
    let settled = false
    let terminate:(() => void)|null = null

    terminate = min.deflate(fixture, (err, buf) => {
        error = err || null
        if (buf) {
            terminate = min.inflate(buf, (err2, buf2) => {
                error = err2 || null
                result = buf2 || null
                settled = true
            })
        } else {
            settled = true
        }
    })

    // Wait for async operations to complete with timeout
    // `settled` is assigned only from the worker callbacks above, which
    // ESLint cannot see, so testing it in the loop header trips
    // no-unmodified-loop-condition. Breaking on it instead is the same
    // check in the same place.
    const start = Date.now()
    while (Date.now() - start < 30000) {
        if (settled) break
        await new Promise(resolve => setTimeout(resolve, 10))
    }

    if (!settled && terminate) terminate()

    t.ok(settled, 'async roundtrip completed')
    t.ok(error === null, `no error: ${error?.message || 'ok'}`)
    t.ok(
        result && result.length === fixture.length &&
      result.every((v, i) => v === fixture[i]),
        'async roundtrip preserves 600000 byte data through worker'
    )
})
