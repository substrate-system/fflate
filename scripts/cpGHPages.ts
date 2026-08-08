import { simpleGit } from 'simple-git'
import { resolve, join } from 'path'
import { cpSync, readdirSync, rmSync, statSync, unlinkSync } from 'fs'

const baseDir = resolve(import.meta.dirname, '..')
const to = (...paths:Array<string>) => join(baseDir, ...paths)
const git = simpleGit()

const branch = (await git.revparse(['--abbrev-ref', 'HEAD'])).trim()
const log = await git.log({ from:'HEAD~1', to:'HEAD' })
const hash = log.latest!.hash.slice(0, 7)

try {
    await git.checkout('gh-pages')

    for (const f of readdirSync(to('.'))) {
        // statSync needs the repository relative path, not the bare entry.
        if (statSync(to(f)).isFile()) unlinkSync(to(f))
    }

    // Remove stale assets before copying new ones, so nested files from
    // past deploys do not accumulate forever.
    rmSync(to('assets'), { recursive:true, force:true })

    // Vite emits nested output, so assets/ is a directory. copyFileSync
    // throws on a directory; cpSync copies the tree. The source lives
    // inside the destination, which node allows, and the walk only visits
    // entries under public/, so this does not recurse into itself.
    const added = readdirSync(to('public'))
    cpSync(to('public'), to('.'), { recursive:true })

    // --update stages every tracked modification and deletion, including
    // the stale files under assets/. It takes no pathspec, so an
    // untracked stray at the root cannot abort the add. The second call
    // stages the new tree; every name in `added` exists after the copy.
    // Neither stages ignored directories, which matters because the
    // branch switch left the working tree without a .gitignore.
    await git.add(['--update'])
    await git.add(added)
    await git.commit('Build demo from ' + hash)
} finally {
    // Always return to the branch we started on, not an assumed master.
    await git.checkout(branch)
}
