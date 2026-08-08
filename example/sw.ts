/// <reference lib="webworker" />

// Required, and not removable. This file has no other import or
// export, so without it TS treats the file as a global script rather
// than a module. In a global script `declare const self` collides with
// lib.dom's `declare var self`, giving TS2451 "Cannot redeclare
// block-scoped variable 'self'" followed by TS2339 on __WB_MANIFEST.
export {}

declare const self:ServiceWorkerGlobalScope&{
    __WB_MANIFEST:Array<{ url:string; revision:string|null }>
}

const manifest = self.__WB_MANIFEST

// vite-plugin-pwa injects per asset revisions but no build version.
// Derive a stable cache name from the revisions so it changes when,
// and only when, the precached set changes.
const precacheVersion = 'fflate-' + manifest
    .map(e => e.revision ?? e.url)
    .join('|')
    .split('')
    .reduce((h, c) => (((h << 5) - h) + c.charCodeAt(0)) | 0, 0)
    .toString(36)

const precacheFiles = manifest.map(e => e.url)

const ch = () => caches.open(precacheVersion)

self.addEventListener('install', ev => {
    // Do not finish installing until every file in the app has been cached
    ev.waitUntil(
        ch().then(
            cache => cache.addAll(precacheFiles)
        )
    )
})

self.addEventListener('activate', ev => {
    ev.waitUntil(
        caches.keys().then(keys => Promise.all(
            keys.filter(k => k !== precacheVersion).map(
                k => caches.delete(k)
            )
        )).then(() => self.clients.claim())
    )
})

self.addEventListener('fetch', ev => {
    ev.respondWith(
        caches.match(ev.request).then(resp => resp || ch().then(c =>
            fetch(ev.request).then(res => c.put(ev.request, res.clone()).then(
                () => res
            ))
        ))
    )
})
