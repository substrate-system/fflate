**What can't you do right now?**

When extracting ZIP archives with fflate, there's no reliable way to distinguish between directories and empty files. Both appear as entries with zero-length data, and the only available heuristic is checking if the filename ends with `/`, which is unreliable since:

- Not all ZIP creation tools add trailing slashes to directory entries
- Some tools create directory entries without the trailing slash convention
- Cross-platform archives may have inconsistent directory representation
- The ZIP specification doesn't require the trailing slash convention

This makes it impossible to properly reconstruct directory structures or handle empty files correctly when extracting archives.

**An optimal solution**

Expose ZIP entry metadata, particularly the external file attributes that contain the directory flag:

```javascript
// Current API only provides filename and data
unzip(zipData, (err, files) => {
  for (const [filename, data] of Object.entries(files)) {
    // Can only guess: filename.endsWith('/') && data.length === 0
  }
});

// Proposed API with metadata
unzip(zipData, (err, files, metadata) => {
  for (const [filename, data] of Object.entries(files)) {
    const meta = metadata[filename];
    if (meta.isDirectory) {
      console.log(`Directory: ${filename}`);
    } else if (data.length === 0) {
      console.log(`Empty file: ${filename}`);
    } else {
      console.log(`File: ${filename}`);
    }
  }
});
```

**(How) is this done by other libraries?**

**JSZip** provides a `dir` property on file objects:
```javascript
zip.forEach((relativePath, file) => {
  if (file.dir) {
    console.log("Directory:", relativePath);
  }
});
```

**node-stream-zip** exposes `isDirectory` directly:
```javascript
const entries = await zip.entries();
for (const entry of Object.values(entries)) {
  if (entry.isDirectory) {
    console.log("Directory:", entry.name);
  }
}
```

**yauzl** (Node.js) provides access to external file attributes:
```javascript
entry.isDirectory = (entry.externalFileAttributes & 0x40000000) !== 0;
```

This feature would enable proper directory handling, which is essential for applications that need to recreate file system structures from ZIP archives or distinguish between intentionally empty files and directory placeholders.