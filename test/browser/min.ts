// The minified bundle ships without declarations, and generating them
// would be meaningless: its public API is identical to the plain
// bundle's by construction.
//
// This file only ISOLATES the untyped import so the cast lives in one
// place. It asserts nothing -- `minBundle` is `any`, so the cast is
// unchecked at compile time and absent at runtime. The real proof that
// the two APIs agree is behavioural: index.ts runs every suite against
// both bundles, and test/3-node-min.ts compares the export-name sets
// on the node side.
// @ts-expect-error no declarations are emitted for the minified bundle
import * as minBundle from '../../dist/browser/index.min.js'
import type { Fflate } from './util.js'

export const min = minBundle as Fflate
