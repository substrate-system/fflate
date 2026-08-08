`FlateError` is exported as an interface rather than a class, which means that
it's not possible to detect it via `instanceof FlateError`.

```js
import { unzipSync, FlateError } from 'fflate'

try {
  unzipSync(...)
} catch (error) {
  if (error instanceof FlateError) {
    ... report the error.code ...
  } else {
    throw error
  }
}
```

>  Exception during run: SyntaxError: The requested module 'fflate' does not
provide an export named 'FlateError'

Hence, it's not possible to determine whether it threw because of the `.zip`
file being "corrupt" or because of an unrelated issue such as
"out of memory", etc.
