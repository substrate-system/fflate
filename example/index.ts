import { zip, unzip, createZippable } from '../src/index.js'
import { humanFilesize } from '@substrate-system/util/filesize'
import { equals } from 'uint8arrays'
import Debug from '@substrate-system/debug'
const debug = Debug()

let pendingDir:FileList
const button = document.querySelector('button[type="submit"')

document.getElementById('dir')?.addEventListener('change', ev => {
    debug('change', ev)
    pendingDir = (ev.target as HTMLInputElement).files!
    debug('pending dir', pendingDir)
    if (pendingDir) {
        button?.removeAttribute('disabled')
    } else {
        button?.setAttribute('disabled', '')
    }
})

document.querySelector('form.example')?.addEventListener('submit', async ev => {
    ev.preventDefault()
    debug('submit', ev)
    const form = (ev.target as HTMLFormElement)
    const files = (form.elements['dir'] as HTMLInputElement).files!

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

            debug('unzipped is equal to pre-zip ???', equals(
                data['abc/test.txt'],
                zippable['abc/test.txt'] as Uint8Array
            ))
        })
    }
})
