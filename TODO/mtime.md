<!-- This template is just a suggestion, feel free to ignore or delete it -->

**What can't you do right now?**
`mtime` of a zip entry is not exposed when unpacking.

**An optimal solution**
* Add `mtime` property to the `filter` function in `unzip` options
* Add `mtime` property to [UnzipFile](https://github.com/101arrowz/fflate/blob/master/docs/interfaces/UnzipFile.md)

Right now I'm patching the library in webpack build to expose time as a 32-bit number `b4(data, o+12)`, converted later on demand (although the calc is trivial so you might just do it by default):
```js
new Date(
  /*Y*/((time >> 25) & 0x7f) + 1980,
  /*M*/(time >> 21) & 0x0f - 1,
  /*D*/(time >> 16) & 0x1f,
  /*h*/(time >> 11) & 0x1f,
  /*m*/(time >> 5) & 0x3f,
  /*s*/(time & 0x1f) * 2);
```

**(How) is this done by other libraries?**
Previously I've used `@zip.js` which exposes it as `lastModDate` that is a JS `Date`.
