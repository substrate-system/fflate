// Runs every suite twice: once against the plain bundle and once
// against the minified one.
//
// The minified pass is the point of this file. The async APIs build
// worker source by parsing identifier names out of array literal TEXT
// via Function.prototype.toString, so a minifier that inlines a
// constant into one of those literals breaks name recovery silently:
// the sync APIs keep working while the async ones corrupt output.
// Nothing else in this repository would catch that.
import * as plainBundle from '../../dist/browser/index.js'
import { min } from './min.js'
import type { Fflate } from './util.js'
import { asyncSuite } from './async.js'
import { streamSuite } from './streams.js'
import { zipSuite } from './zip.js'
import { harnessSuite } from './harness.js'

const plain = plainBundle as unknown as Fflate

// Runs once, not per bundle: it checks the runner's own configuration,
// which is what decides whether the minified pass below runs at all.
harnessSuite()

for (const [label, f] of [
    ['plain', plain],
    ['min', min]
] as const) {
    asyncSuite(f, label)
    streamSuite(f, label)
    zipSuite(f, label)
}
