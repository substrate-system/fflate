import * as esbuild from 'esbuild'
import { minify, type MinifyOptions } from 'terser'
import { execFileSync } from 'child_process'
import { rmSync, mkdirSync, writeFileSync, readFileSync } from 'fs'
import { resolve, join } from 'path'

const root = resolve(import.meta.dirname, '..')
function p (...parts:Array<string>):string {
    return join(root, ...parts)
}

// Step 1: remove any output from a previous run.
rmSync(p('dist'), { recursive:true, force:true })
mkdirSync(p('dist/browser'), { recursive:true })
mkdirSync(p('dist/node'), { recursive:true })
mkdirSync(p('dist/umd'), { recursive:true })

console.log('cleaned dist/')

// Proven against wcln by the existing UMD build. Do not add compress
// options without re-running the Phase 6 browser suite against the
// minified bundle: any option that substitutes a constant's VALUE for
// its NAME inside the bInflt/bDflt array literals breaks name recovery
// silently, leaving sync APIs working while async APIs corrupt output.
function terserOpts (esm:boolean):MinifyOptions {
    return {
        mangle:{ toplevel:true },
        compress:{ passes:5, unsafe:true, pure_getters:true },
        module:esm
    }
}

type Bundle = { code:string; map:string }

// esbuild with write:false returns the bundle and its map as separate
// entries in outputFiles. Confirmed against esbuild 0.25.12.
function collect (out:esbuild.BuildResult):Bundle {
    const js = out.outputFiles!.find(f => !f.path.endsWith('.map'))!
    const map = out.outputFiles!.find(f => f.path.endsWith('.map'))!
    return { code:js.text, map:map.text }
}

async function writeBundle (
    dir:string,
    name:string,
    bundle:Bundle,
    esm:boolean
):Promise<void> {
    writeFileSync(p(dir, name + '.js'), bundle.code)
    writeFileSync(p(dir, name + '.js.map'), bundle.map)

    const min = await minify(bundle.code, {
        ...terserOpts(esm),
        sourceMap:{
            content:bundle.map,
            url:name + '.min.js.map'
        }
    })
    if (!min.code) throw new Error('terser produced no output for ' + name)

    writeFileSync(p(dir, name + '.min.js'), min.code)
    // terser types map as string|RawSourceMap, but it returns a string
    // whenever sourceMap.content is supplied, which it always is here.
    writeFileSync(p(dir, name + '.min.js.map'), min.map as string)

    console.log('wrote ' + dir + '/' + name + '{.js,.min.js} and maps')
}

// Step 2: declarations. Invoke the local tsc through node rather than
// npx: npx resolves to npx.cmd on Windows, which execFileSync will not
// find without a shell.
execFileSync(
    process.execPath,
    [
        p('node_modules/typescript/bin/tsc'),
        '--emitDeclarationOnly',
        '--project',
        'tsconfig.build.json'
    ],
    { cwd:root, stdio:'inherit' }
)

const dts = readFileSync(p('dist/index.d.ts'), 'utf8')

// The bundled builds expose no separate worker entry point, so only
// index.d.ts is published. Declaration maps are dropped with it: they
// would point at src/ paths that are not shipped.
for (const f of [
    'index.d.ts', 'index.d.ts.map',
    'worker.d.ts', 'worker.d.ts.map',
    'node-worker.d.ts', 'node-worker.d.ts.map'
]) {
    rmSync(p('dist', f), { force:true })
}

writeFileSync(p('dist/browser/index.d.ts'), dts)
writeFileSync(p('dist/node/index.d.ts'), dts)

console.log('wrote declarations')

// Step 3: browser bundle.
// esbuild's --alias only rewrites bare specifiers, so the relative
// ./node-worker to ./worker swap needs a resolve plugin.
const workerSwap:esbuild.Plugin = {
    name:'browser-worker-swap',
    setup (build) {
        build.onResolve({ filter:/^\.\/node-worker$/ }, args => ({
            path:resolve(args.resolveDir, 'worker.ts')
        }))
    }
}

const browser = collect(await esbuild.build({
    entryPoints:[p('src/index.ts')],
    outfile:p('dist/browser/index.js'),
    bundle:true,
    write:false,
    format:'esm',
    platform:'browser',
    target:'es2022',
    // keepNames must stay off. It wraps every function as
    // __name(fn, 'name'), and the async APIs stringify functions and eval
    // them inside a worker where that helper is not in scope, giving
    // "ReferenceError: __name is not defined". Name recovery for the
    // worker payload is wcln's job, not esbuild's.
    keepNames:false,
    sourcemap:true,
    plugins:[workerSwap]
}))

await writeBundle('dist/browser', 'index', browser, true)

// Step 4: node bundle.
const node = collect(await esbuild.build({
    entryPoints:[p('src/index.ts')],
    outfile:p('dist/node/index.js'),
    bundle:true,
    write:false,
    format:'esm',
    platform:'node',
    target:'es2022',
    // Off for the same reason as the browser bundle above.
    keepNames:false,
    sourcemap:true,
    external:['node:worker_threads']
}))

await writeBundle('dist/node', 'index', node, true)

// Step 5: UMD. Bundles the browser variant, so the swap plugin applies.
const umd = await esbuild.build({
    entryPoints:[p('src/index.ts')],
    outfile:p('dist/umd/fflate.js'),
    bundle:true,
    write:false,
    format:'iife',
    globalName:'fflate',
    platform:'browser',
    target:'es2022',
    // Off for the same reason as the browser bundle: async APIs
    // stringify functions and eval them in a worker context where
    // __name is not in scope.
    keepNames:false,
    sourcemap:false,
    plugins:[workerSwap]
})

const iife = umd.outputFiles!.find(f => !f.path.endsWith('.map'))!.text

const wrapped = [
    '(function(root, factory){',
    "  if (typeof module === 'object' && typeof exports === 'object')",
    '    module.exports = factory();',
    "  else if (typeof define === 'function' && define.amd)",
    '    define([], factory);',
    '  else root.fflate = factory();',
    "})(typeof self !== 'undefined' ? self : this, function(){",
    iife,
    '  return fflate;',
    '});'
].join('\n')

// Step 6: minify. UMD is CommonJS shaped, so module:false here.
const umdMin = await minify(wrapped, terserOpts(false))
if (!umdMin.code) throw new Error('terser produced no output for umd')

writeFileSync(p('dist/umd/fflate.js'), umdMin.code)

// The root package.json is "type": "module", so node would parse this
// .js as ESM, where the wrapper's `this` is undefined and loading throws
// "Cannot set properties of undefined". Scoping the directory back to
// commonjs makes the artifact loadable by both require() and import()
// while keeping the filename that unpkg points at.
writeFileSync(p('dist/umd/package.json'), '{ "type": "commonjs" }\n')

console.log('wrote dist/umd/fflate.js')
