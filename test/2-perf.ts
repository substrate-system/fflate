import type { TestHandler } from './util.js'
import { testSuites, workers, bClone } from './util.js'
import { writeFileSync } from 'fs'
import { join } from 'path'
import { fileURLToPath } from 'url'

const here = fileURLToPath(new URL('.', import.meta.url))

const preprocessors = {
    inflate:workers.zlib.deflate,
    gunzip:workers.zlib.gzip,
    unzlib:workers.zlib.zlib
}

const cache:Record<string, Record<string, Buffer>> = {
    deflate:{},
    inflate:{},
    gzip:{},
    gunzip:{},
    zlib:{},
    unzlib:{}
}

const flattenedWorkers:Record<string, TestHandler> = {}
for (const k in workers) {
    for (const l in workers[k as keyof typeof workers]) {
        if (l === 'zip' || l === 'unzip') continue
        flattenedWorkers[k + '.' + l] = async (file, name, resetTimer) => {
            const fileClone = bClone(file)
            let buf = fileClone
            if (preprocessors[l as keyof typeof preprocessors]) {
                buf = bClone(cache[l][name] ||= Buffer.from(
                    await preprocessors[l as keyof typeof preprocessors](buf, [buf.buffer])
                ))
                resetTimer()
            }
            // undefined rather than null: pako 3 dereferences the options
            // argument directly, so an explicit null throws.
            const opt2 = preprocessors[l as keyof typeof preprocessors]
                ? k === 'tinyInflate'
                    ? new Uint8Array(file.length)
                    : undefined
                : { level:1 }
            await workers[k as keyof typeof workers][l as 'inflate']([buf, opt2], opt2 instanceof Uint8Array
                ? [buf.buffer, opt2.buffer]
                : [buf.buffer])
        }
    }
}

testSuites(flattenedWorkers).then(perf => {
    writeFileSync(join(here, 'results', 'timings.json'), JSON.stringify(perf, null, 2))
})
