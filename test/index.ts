// Entry point for the node test suite.
//
// tapzero has no CLI runner: importing a suite registers its tests, and
// the runner flushes them in registration order once this module
// finishes evaluating. The numeric filename prefixes therefore set the
// run order.
//
// ZIP, stream, and async coverage lives in test/browser/, where it runs
// against the real browser build including the minified bundle.
import './0-valid.js'
import './1-size.js'
import './2-perf.js'
import './3-node-min.js'
import './4-error-class.js'
import './5-async-error.js'
import './6-unzip-attrs.js'
import './7-unzip-rawname.js'
import './8-stream-chunks.js'
import './9-flush-sync.js'
import './10-stream-window.js'
