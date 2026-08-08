<!-- This template is just a suggestion, feel free to ignore or delete it -->
**How to reproduce**

```js
import { inflateRawSync } from 'node:zlib'
import { Deflate } from 'fflate'

function compressRawWithFlush(input, syncFlush) {
  const chunks = []
  const deflate = new Deflate((chunk) => chunks.push(chunk))

  deflate.push(new TextEncoder().encode(input))
  deflate.flush(syncFlush)
  deflate.push(new Uint8Array(0), true)

  return Buffer.concat(chunks)
}

const input = 'a'.repeat(70000)

const withSyncFlush = compressRawWithFlush(input, true)
const withoutSyncFlush = compressRawWithFlush(input, false)

try {
  inflateRawSync(withSyncFlush)
  console.log('sync flush: OK')
} catch (error) {
  console.log('sync flush: FAIL ->', error.message)
}

try {
  const output = inflateRawSync(withoutSyncFlush)
  console.log('non-sync flush: OK ->', output.length, 'bytes')
} catch (error) {
  console.log('non-sync flush: FAIL ->', error.message)
}
```

```
sync flush: FAIL -> invalid stored block lengths
non-sync flush: OK -> 70000 bytes
```

**The problem**

Compressed data using sync flush fails to inflate. This is dependent on input: small input don't show this issue.
