import assert from 'assert'
import { zip, unzip, createZippable } from '../src/index.js'
import { humanFilesize } from '@substrate-system/util/filesize'
import { equals } from 'uint8arrays'
import Debug from '@substrate-system/debug'
const debug = Debug()

document.getElementById('file')?.addEventListener('change', ev => {
    debug('change', ev)
})

document.querySelector('form.example')?.addEventListener('submit', async ev => {
    ev.preventDefault()
    debug('submit', ev)
    const form = (ev.target as HTMLFormElement)
    const files = (form['file'] as HTMLInputElement).files!
    
    // const zippable = await Array.from(files).reduce(async (_acc, file) => {
    //     const acc = await _acc
    //     acc[file.webkitRelativePath] = new Uint8Array(await file.arrayBuffer())
    //     return acc
    // }, Promise.resolve({}) as Promise<Record<string, Uint8Array>>)

    const zippable = await createZippable(files)

    debug('zippable', zippable)

    const preZipSize = Object.keys(zippable).reduce((total, key) => {
        return total + (zippable[key] as Uint8Array).length
    }, 0)
    debug('size, before zipping...', humanFilesize(preZipSize))

    zip(zippable, {
        level: 6
    }, (err, data) => {
        if (err) throw err
        debug('size, after zipping...', humanFilesize(data.length))

        // now unzip
        exampleUnzip(data)
    })

    function exampleUnzip (zippedData:Uint8Array) {
        unzip(zippedData, (err, data) => {
            if (err) throw err
            debug('the unzipped data...', data)

            debug('unzipped equals zipped???', equals(
                data['abc/test.txt'],
                zippable['abc/test.txt'] as Uint8Array
            ))
        })
    }
})
