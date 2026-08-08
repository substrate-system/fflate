## Summary

Between 0.8.2 and 0.8.3, the streaming `Gunzip` callback's emission pattern changed. For small inputs that 0.8.2 produced in a single callback, 0.8.3 now emits the data chunk followed by **two empty chunks** (one with `final: false`, one with `final: true`).

The documented contract for `FlateStreamHandler` doesn't forbid this, so the change is arguably within spec. But it's a behavior change that breaks downstream code that captured "the chunk" rather than concatenating — at least [`@arethetypeswrong/core`](https://github.com/arethetypeswrong/arethetypeswrong.github.io/blob/main/packages/core/src/createPackage.ts) does this, and `attw -P` / `attw --from-npm` is now broken for every package as a result.

I suspect this is a side-effect of the new "Support sync flushes (`Z_SYNC_FLUSH` in zlib)" feature in the 0.8.3 changelog, but I wanted to confirm and raise a couple of questions.

## Repro

Two clean installs, same 11-byte input:

```sh
mkdir /tmp/f82 /tmp/f83
cd /tmp/f82 && npm init -y && npm i fflate@0.8.2
cd /tmp/f83 && npm init -y && npm i fflate@0.8.3
node -e '
  const f82 = require("/tmp/f82/node_modules/fflate");
  const f83 = require("/tmp/f83/node_modules/fflate");
  const gz = f82.gzipSync(new TextEncoder().encode("hello world"));
  for (const [name, F] of [["0.8.2", f82], ["0.8.3", f83]]) {
    const cbs = [];
    new F.Gunzip((c, final) => cbs.push({ len: c.length, final })).push(gz, true);
    console.log(name + ":", JSON.stringify(cbs));
  }
'
```

Output:

```
0.8.2: [{"len":11,"final":true}]
0.8.3: [{"len":11,"final":false},{"len":0,"final":false},{"len":0,"final":true}]
```

## Questions

1. **Is this intentional?** It looks like a consequence of the sync-flush work in 0.8.3. If so, that's reasonable, but it would be useful to note in the changelog/docs that the streaming callback may now fire with `final: false` for inputs that previously emitted a single final block.

2. **Should empty chunks be suppressed?** Emitting `{ len: 0, final: false }` followed by `{ len: 0, final: true }` is the most surprising part for consumers — both chunks have no data. If the intent is to signal "flush boundary" and "end of stream," the empty `final: false` chunk in particular seems like it could be skipped without losing information.

3. **If this is intended and considered non-breaking,** would you be open to a documentation note on `FlateStreamHandler` clarifying the new emission pattern? Downstream code that relied on "single callback for small inputs" was making an undocumented assumption, but the assumption was widely safe under 0.8.2.

## Downstream impact

`@arethetypeswrong/core@0.18.2` does this (paraphrased):

```js
let unzipped;
new Gunzip((chunk) => (unzipped = chunk)).push(tarball, /*final*/ true);
const data = untar(unzipped);
```

Under 0.8.2, `unzipped` is the full inflated tarball. Under 0.8.3, it's the empty final-flush chunk, so `untar([])` returns `[]` and the next line throws. Result: `attw` is broken on every install today because attw declares `fflate: "^0.8.2"`. Workaround on the attw side is straightforward (accumulate chunks), and I've filed that separately, but I wanted to surface the fflate-side change since it may affect other consumers similarly. See https://github.com/arethetypeswrong/arethetypeswrong.github.io/issues/262

## Environment

- fflate 0.8.2 (works) vs 0.8.3 (changed)
- Node 26.2.0, macOS 25.5.0
