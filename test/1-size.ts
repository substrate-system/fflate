import { testSuites, workers, bClone } from './util.js'
import { writeFileSync } from 'fs'
import { join } from 'path'
import { fileURLToPath } from 'url'
import { performance } from 'perf_hooks'

const here = fileURLToPath(new URL('.', import.meta.url))

const sizePerf:Record<string, Record<string, [number, number]>> = {}

testSuites({
    async main (file, name, _resetTimer, t) {
        sizePerf[name] = {}
        for (const lib of (['fflate', 'pako', 'uzip', 'zlib'] as const)) {
            const clone = bClone(file)
            const ts = performance.now()
            sizePerf[name][lib] = [(await workers[lib].deflate([clone, { level:9 }], [clone.buffer])).length, performance.now() - ts]
        }
        for (const lib of ['pako', 'uzip', 'zlib']) {
            // Less than 5% larger
            t.ok(
                ((sizePerf[name].fflate[0] - sizePerf[name][lib][0]) / sizePerf[name][lib][0]) < 0.05,
                'fflate should be within 5% of ' + lib + ' on ' + name
            )
        }
    }
}).then(() => {
    writeFileSync(join(here, 'results', 'longTimings.json'), JSON.stringify(sizePerf, null, 2))
})
