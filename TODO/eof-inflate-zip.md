

<!-- If possible, upload files or instructions to reproduce the bug -->
<!-- You can send them to me privately at arjunbarrett@gmail.com if they contain confidential information -->

**The problem**

I have a zip file that I can decompress with other decompressors (ubuntu's `unzip`, yazul, python3) but when using fflate the zip is unable to be decompressed.

```
Error: unexpected EOF
    at inflt ([worker eval]:233:25)
    at Inflate.c ([worker eval]:289:18)
    at Inflate.push ([worker eval]:294:29)
    at [worker eval]:263:40
    at MessagePort.<anonymous> ([worker eval]:298:86)
```

<!-- This template is just a suggestion, feel free to ignore or delete it -->
**How to reproduce**

Using a simple script

```javascript
import { AsyncUnzipInflate, Unzip } from 'fflate';

const bytes = readFileSync('./address.zip'); 
const unzip = new Unzip((file) => {
  console.log(file.name);
  file.ondata = (err, data, final) => {
    if (err) throw err;
    console.log('got', file.name, data.length, final);
  };
  file.start();
});

unzip.register(AsyncUnzipInflate);
unzip.push(bytes);
```

with my sample zip file: https://public.chard.com/fflate/address.zip