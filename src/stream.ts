import { Zip, ZipPassThrough } from './index.js'

const wrapError = err => {
    if (err instanceof Error) {
        return err
    } else if (typeof err === 'string') {
        return new Error(`String error from fflate: ${err}`)
    } else {
        return new Error(`Unexpected error type from fflate: ${err}`)
    }
}

export function makeZipStream (files:File[]) {
    const highWaterMark = 1_000_000

    const zip = new Zip()

    // Input state
    const compressingQueue:number[] = []  // Input size of each chunk inside compressor
    let numCompressingBytes = 0  // Total size inside compressor

    let waitingInput:(null|{ resolve, reject }) = null
    const waitForSpace = async () => {
        while (
            // eslint-disable-next-line no-unmodified-loop-condition
            !outputWantsData ||
            // eslint-disable-next-line no-unmodified-loop-condition
            numCompressingBytes >= highWaterMark
        ) {
            try {
                await new Promise((resolve, reject) => {
                    waitingInput = {
                        resolve,
                        reject
                    }
                })
            } finally {
                waitingInput = null
            }
        }
    }

    // Output state
    let aborted = false
    let outputWantsData = false

    // Updated for each file
    let fileReader:null|ReadableStreamDefaultReader<Uint8Array> = null
    let readableController:(null|ReadableStreamDefaultController<ArrayBufferView>) = null

    const readable = new globalThis.ReadableStream(
        {
            start (controller) {
                readableController = controller
                zip.ondata = (err, chunk, final) => {
                    if (err) {
                        aborted = true
                        fileReader?.cancel(wrapError(err))
                        return
                    }

                    controller.enqueue(chunk)
                    if (final) {
                        controller.close()
                        return
                    }

                    if (controller.desiredSize! <= 0) {
                        outputWantsData = false
                    }
                }
            },
            pull () {
                outputWantsData = true
                waitingInput?.resolve()
            },
            cancel (reason) {
                aborted = true
                fileReader?.cancel(reason)
            }
        },
        new globalThis.ByteLengthQueuingStrategy({
            highWaterMark
        })
    )

  ;(async () => {
        for (const file of files) {
            if (aborted) {
                return
            }

            const { stream } = file
            const path = file.webkitRelativePath
            const zipItem = new ZipPassThrough(path)
            zip.add(zipItem)

            const originalOnData = zipItem.ondata.bind(zipItem)
            zipItem.ondata = (err, chunk, isLast) => {
                originalOnData(err, chunk, isLast)

                numCompressingBytes -= Number(compressingQueue.shift())
                if (err) {
                    waitingInput?.reject(wrapError(err))
                } else {
                    waitingInput?.resolve()
                }
            }

            try {
                fileReader = stream().getReader()
                let fileDone = false
                while (!fileDone) {
                    const { done, value } = await fileReader.read()
                    const chunk = done ? new Uint8Array(0) : value
                    fileDone = done

                    await waitForSpace()

                    compressingQueue.push(chunk.byteLength)
                    numCompressingBytes += chunk.byteLength
                    zipItem.push(chunk, done)
                }
            } catch (err) {
                readableController!.error(err)
            }
        }

        zip.end()
    })()

    return readable
}
