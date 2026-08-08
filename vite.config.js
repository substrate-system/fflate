import { defineConfig } from 'vite'
import preact from '@preact/preset-vite'
import { VitePWA } from 'vite-plugin-pwa'
import { resolve } from 'path'

const root = import.meta.dirname

export default defineConfig({
    root: 'example',
    publicDir: '_public',

    optimizeDeps: {
        // src/index.ts has @ts-ignore on method signatures without
        // implementation, assigned at runtime by the worker wrapper.
        // rolldown's scanner does not honor @ts-ignore, so disable
        // pre-bundling for fflate to avoid parse errors in dev.
        exclude: ['fflate']
    },

    resolve: {
        // Array form, because the worker swap below matches on a regex.
        alias: [
            // preset-vite aliases react, react-dom, react/jsx-runtime
            // and react-dom/test-utils, but not the react-dom/client
            // subpath that example/index.tsx imports createRoot from.
            // The target is preact/compat/client, not preact/compat:
            // createRoot lives in the client entry only.
            { find: 'react-dom/client', replacement: 'preact/compat/client' },

            // The example imports the library as '../../..', which would
            // resolve through the root package.json and couple this
            // build to the library build. Resolve it to source instead.
            { find: 'fflate', replacement: resolve(root, 'src/index.ts') },

            // Because the line above pulls in src/index.ts, this build
            // sees its `import wk from './node-worker'`. From Phase 2
            // onward that module statically imports node:worker_threads,
            // which must never reach a browser bundle. Apply the same
            // swap the library build performs.
            {
                find: /^\.\/node-worker$/,
                replacement: resolve(root, 'src/worker.ts')
            }
        ]
    },

    plugins: [
        preact(),
        VitePWA({
            strategies: 'injectManifest',
            srcDir: '.',
            filename: 'sw.ts',
            injectRegister: false,
            manifest: false,
            injectManifest: {
                injectionPoint: 'self.__WB_MANIFEST'
            }
        })
    ],

    build: {
        // outDir is outside root, so Vite refuses to clear it unless
        // emptyOutDir is set explicitly.
        outDir: resolve(root, 'public'),
        emptyOutDir: true
    }
})
