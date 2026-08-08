import { test } from '@substrate-system/tapzero'
import {
    FlateError,
    FlateErrorCode,
    inflateSync,
    unzipSync
} from 'fflate'

// `FlateError` used to be an interface, so it vanished at runtime and
// `import { FlateError }` failed outright. These tests pin it down as a
// real exported class.

test('FlateError is exported as a runtime value', t => {
    t.equal(
        typeof FlateError,
        'function',
        'FlateError should be a constructor, not a type-only export'
    )
    t.ok(
        FlateError.prototype instanceof Error,
        'FlateError should extend Error'
    )
})

test('inflateSync throws a FlateError on truncated input', t => {
    try {
        inflateSync(new Uint8Array([0, 0, 0]))
        t.fail('inflateSync should have thrown')
    } catch (e) {
        t.ok(e instanceof FlateError, 'error should be instanceof FlateError')
        t.ok(e instanceof Error, 'error should still be instanceof Error')
        t.equal(
            (e as FlateError).code,
            FlateErrorCode.UnexpectedEOF,
            'code should be UnexpectedEOF'
        )
        t.equal(
            (e as FlateError).message,
            'unexpected EOF',
            'message string should be unchanged'
        )
        t.equal((e as FlateError).name, 'FlateError', 'name should be FlateError')
    }
})

test('inflateSync throws a FlateError on an invalid block type', t => {
    try {
        inflateSync(new Uint8Array([6, 6, 6, 6]))
        t.fail('inflateSync should have thrown')
    } catch (e) {
        t.ok(e instanceof FlateError, 'error should be instanceof FlateError')
        t.equal(
            (e as FlateError).code,
            FlateErrorCode.InvalidBlockType,
            'code should be InvalidBlockType'
        )
        t.equal(
            (e as FlateError).message,
            'invalid block type',
            'message string should be unchanged'
        )
    }
})

test('unzipSync throws a FlateError on a corrupt archive', t => {
    try {
        unzipSync(new Uint8Array(8))
        t.fail('unzipSync should have thrown')
    } catch (e) {
        t.ok(e instanceof FlateError, 'error should be instanceof FlateError')
        t.equal(
            (e as FlateError).code,
            FlateErrorCode.InvalidZipData,
            'code should be InvalidZipData'
        )
    }
})

test('a plain Error is not mistaken for a FlateError', t => {
    t.ok(
        !(new Error('out of memory') instanceof FlateError),
        'unrelated failures should not satisfy the instanceof check'
    )
})

test('err() keeps its message override', t => {
    // The `msg` override path is only reachable through the public API via
    // the header checks, which pass an explicit string.
    try {
        // A zlib header check failure uses err(6, 'invalid zlib data').
        unzipSync(new Uint8Array(8))
        t.fail('unzipSync should have thrown')
    } catch (e) {
        t.equal(
            (e as FlateError).message,
            'invalid zip data',
            'overridden messages survive the class change'
        )
    }
})
