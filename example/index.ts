import { zip } from '../src/index.js'
import Debug from '@substrate-system/debug'
const debug = Debug()

document.getElementById('file')?.addEventListener('change', ev => {
    debug('change', ev)
})

document.querySelector('form.example')?.addEventListener('submit', ev => {
    ev.preventDefault()
    debug('submit', ev)
    const form = (ev.target as HTMLFormElement)
    const files = (form['file'] as HTMLInputElement).files!
    debug('the files', files)
    
    const zippable = Array.from(files).reduce((acc, file) => {
        acc[file.webkitRelativePath] = file
        return acc
    }, {})

    debug('zippable', zippable)

    zip(zippable, {
        level: 6
    }, (err, data) => {
        debug('zipped the things...', err)
        debug('the data...', data)
    })
})
