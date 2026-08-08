import { test } from '@substrate-system/tapzero'
import pkg from '../../package.json'
import {
    DEFAULT_TIMEOUT_MS, TIMEOUT_MARGIN_MS, autoFinishWindow
} from './util.js'

// Guards the one invariant that, when broken, makes this whole suite
// lie rather than fail.
//
// tapout ends a run after a window of console silence and resets that
// timer on every line printed. The suites are silent while polling, so
// a stalled test goes quiet, auto-finish fires, and every later test is
// dropped -- including the entire minified pass, which is the only
// reason this directory exists. Measured with a deliberately stalled
// test at tapout's default timeout: zero `# min >` headers, no plan
// line, and the run reported as a PASS.
//
// The fix is a pair of numbers in two different files: the timeout flag
// in package.json's test:browser script, and withTimeout's default in
// util.ts. The rejection has to land before auto-finish, because the
// `not ok` it prints is itself what resets the window and lets the run
// continue. Editing either value alone silently restores the
// truncation, so assert the relationship rather than either number.
export function harnessSuite () {
    test('harness > timeout budget leaves room to report a stall', t => {
        const script = (pkg as { scripts:Record<string, string> })
            .scripts['test:browser']

        // tapout accepts `-t` as a full alias for `--timeout`
        // (dist/cli.js:52), and takes the LAST occurrence when a flag is
        // repeated, so match globally and keep the final capture. Reading
        // the first would report green on
        // `--timeout 30000 --timeout 4000`, which is the one input where a
        // naive guard passes a genuinely broken config.
        const flag = /(?:--timeout|-t)[= ]\s*(\d+)/g
        let last:RegExpExecArray|null = null
        for (
            let m = flag.exec(script); m !== null; m = flag.exec(script)
        ) last = m

        t.ok(last !== null, last ?
      `test:browser sets an explicit timeout of ${last[1]}ms` :
            'test:browser must pass an explicit timeout ' +
      '(--timeout <ms> or -t <ms>); tapout\'s default of 5000 puts ' +
      'the auto-finish window below withTimeout\'s ' +
      `${DEFAULT_TIMEOUT_MS}ms default. Script: ${script}`)

        if (!last) return

        // Not a bare `<`. The auto-finish timer starts at the previous
        // test's output line, while withTimeout's starts later, after the
        // test body's setup -- deflateSync over a 600000-byte fixture in
        // several of these. A default just under the window would satisfy
        // `<` and still lose the race.
        const window = autoFinishWindow(Number(last[1]))
        const fits = DEFAULT_TIMEOUT_MS + TIMEOUT_MARGIN_MS <= window
        t.ok(fits, fits ?
      `withTimeout default ${DEFAULT_TIMEOUT_MS}ms fits inside the ` +
      `${window}ms auto-finish window` :
      `withTimeout default (${DEFAULT_TIMEOUT_MS}ms) plus ` +
      `${TIMEOUT_MARGIN_MS}ms of setup margin must fit inside ` +
      `tapout's auto-finish window (${window}ms from timeout ` +
      `${last[1]}), or a stalled test truncates the run instead of ` +
      'reporting')
    })
}
