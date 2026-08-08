**What can't you do right now?**
It happens that in Russia file names inside zip files are often encoded with cp866. Such filenames currently decoded incorrectly in fflate. The best I can do is

```js
  new TextDecoder('cp866').decode(strToU8(file.name))
```

but it produces correct characters interleaved with some gibberish.

**An optimal solution**
Either provide the raw name in UnzipFile
```
{
    name: string,//as it is decoded now
    rawName: {
        bytes: Uint8Array,
        isUTF8: boolean
    },
    ondata: AsyncFlateStreamHandler,
    ...
}
```

, or make it possible to provide an encoding for entries marked as not utf-8.
```js
unzip = new Unzip();
unzip.setFallbackEncoding('cp866');
```

**(How) is this done by other libraries?**
jszip also fails to decode it correctly.

There is `unzip -O cp866` in Ubuntu starting from some version, and before that version I believe they had a hack that would have used cp866 automatically if it had seen a Russian locale in the OS.
A browser equivalent for that hack would be `navigator.language == 'ru-RU'` if you are willing to use that approach.
