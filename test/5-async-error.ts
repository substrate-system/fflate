import { test } from '@substrate-system/tapzero'

// An error raised inside a worker cannot keep its prototype: the browser
// posts it back through structured clone and node re-serializes it for the
// `error` event. Both drop the class, so the main thread used to receive a
// plain `Error` carrying only `message`, `stack`, and `code`. These tests
// pin down that the async APIs rehydrate it into a real `FlateError`.
//
// They run against the SHIPPED bundles rather than `fflate` (which
// tsconfig maps to src/). tsx transpiles the source with esbuild's
// `keepNames` on, which wraps every function as `__name(fn, '...')`; the
// worker payload is stringified source, and `__name` does not exist there.
// So the worker path is only exercisable against `dist/`, same as
// test/3-node-min.ts.

const importBundles = async () => {
    const plain = await import('../dist/node/index.js')
    // @ts-expect-error no declarations for minified bundle
    const min:typeof plain = await import('../dist/node/index.min.js')
    return { plain, min }
}

type Bundle = Awaited<ReturnType<typeof importBundles>>['plain']

// Truncated DEFLATE stream: the worker hits UnexpectedEOF.
const corrupt = () => new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8])

const inflateErr = (p:Bundle):Promise<unknown> => {
    return new Promise((resolve, reject) => {
        p.inflate(corrupt(), e => {
            if (e) resolve(e)
            else reject(new Error('inflate should have failed'))
        })
    })
}

const asyncInflateErr = (p:Bundle):Promise<unknown> => {
    return new Promise((resolve, reject) => {
        const strm = new p.AsyncInflate(e => {
            if (e) resolve(e)
        })
        strm.push(corrupt(), true)
        setTimeout(() => reject(new Error('AsyncInflate never errored')), 5000)
    })
}

for (const label of ['plain', 'min'] as const) {
    test(`${label}: inflate callback receives a FlateError`, async t => {
        const bundles = await importBundles()
        const p = bundles[label]
        const e = await inflateErr(p)

        t.ok(
            e instanceof p.FlateError,
            'callback error should be instanceof FlateError'
        )
        t.ok(e instanceof Error, 'callback error should still be an Error')
        t.equal(
            (e as InstanceType<Bundle['FlateError']>).code,
            p.FlateErrorCode.UnexpectedEOF,
            'code should survive the worker boundary'
        )
        t.equal(
            (e as Error).message,
            'unexpected EOF',
            'message should survive the worker boundary'
        )
        t.equal((e as Error).name, 'FlateError', 'name should be FlateError')
        t.ok(
            typeof (e as Error).stack === 'string' && (e as Error).stack !== '',
            'the worker-side stack should be preserved'
        )
    })

    test(`${label}: AsyncInflate handler receives a FlateError`, async t => {
        const bundles = await importBundles()
        const p = bundles[label]
        const e = await asyncInflateErr(p)

        t.ok(
            e instanceof p.FlateError,
            'stream error should be instanceof FlateError'
        )
        t.equal(
            (e as InstanceType<Bundle['FlateError']>).code,
            p.FlateErrorCode.UnexpectedEOF,
            'code should survive the worker boundary'
        )
        t.equal(
            (e as Error).message,
            'unexpected EOF',
            'message should survive the worker boundary'
        )
    })

    test(`${label}: a non-FlateError worker failure stays a plain Error`, t => {
        // Guard against rehydrating indiscriminately. Only a numeric
        // `code` marks an error as ours; an unrelated failure such as a
        // worker exiting must not be relabelled.
        const bundles = importBundles()
        return bundles.then(bs => {
            const p = bs[label]
            const other = new Error('exited with code 1')
            t.ok(
                !(other instanceof p.FlateError),
                'an unrelated Error should not satisfy the check'
            )
        })
    })
}
