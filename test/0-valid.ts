import { testSuites, workers, bClone } from './util.js'
import { test } from '@substrate-system/tapzero'

// Name is to ensure that this runs first
// Note that workers are not used here to optimize performance but rather
// to prevent infinite loops from hanging the process.
testSuites({
    async compression (file, _name, _resetTimer, t) {
        const fileClone = bClone(file)
        const cProm = workers.fflate.deflate(fileClone, [fileClone.buffer])
        cProm.timeout(10000)
        const buf = await cProm
        t.ok(
            file.equals(await workers.zlib.inflate(buf, [buf.buffer])),
            'fflate deflate output should inflate back to the original'
        )
    },
    async decompression (file, _name, _resetTimer, t) {
        const fileClone = bClone(file)
        const data = await workers.zlib.deflate(fileClone, [fileClone.buffer])
        const dProm = workers.fflate.inflate(data, [data.buffer])
        dProm.timeout(5000)
        t.ok(
            file.equals(await dProm),
            'fflate should inflate zlib output back to the original'
        )
    }
})

// Test error-marshalling path: worker throws, structured error is
// reconstructed as Error instance with stack trace
test('worker error handling', async t => {
    const badData = new Uint8Array([9, 9, 9, 9, 9, 9, 9, 9])
    const promise = workers.fflate.inflate(badData, [badData.buffer])
    promise.timeout(5000)
    try {
        await promise
        t.ok(false, 'should have thrown')
    } catch (e) {
        t.ok(e instanceof Error, 'error should be Error instance')
        t.ok(
            typeof e.message === 'string' && e.message.length > 0,
            'error should have non-empty message string'
        )
        // Not just "is a string" -- that is true of any Error. The point of
        // the marshalling is that the stack comes from the worker, so it
        // must not have been constructed locally in wc().
        t.ok(
            typeof e.stack === 'string' && !e.stack.includes('test/util.ts'),
            'error stack should come from the worker, not from wc()'
        )
    }
})
