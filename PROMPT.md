study TODO/*.md

# MISSION

Fix all documented issues in `fflate`.

## Important

- Commit frequently.
- Read the `## Codebase Patterns` section at the top of `progress.log`
  before starting.
- Complete exactly ONE task this session, then stop.

# EXECUTION RULES

1. READ: At the start of every session, read `specs/prd.json`, `CLAUDE.md`
   (its `## Invariants` and `## Gotchas` sections are load-bearing --
   several of them fail silently), `progress.log` (especially the
   `## Codebase Patterns` section), and the `TODO/*.md` note behind your
   story if one exists.
2. SCOPE: The loop pins your task for you. Work on the story named in the
   `# THIS SESSION` block appended at the end of this prompt, and only that
   one. If there is no such block (you are being run by hand), pick the
   highest-priority story where `passes: false` and `manual` is not `true`.
   Never touch a story flagged `"manual": true` -- those need a human at a
   real browser and are handled by `./ralph_device.sh`. If the task has
   acceptance criteria, treat them as the definition of done; if it doesn't,
   infer the smallest change that fully satisfies the title and note your
   interpretation in `progress.log`.
3. TEST-FIRST (where applicable): For a feature or bug fix, write failing tests
   that capture the desired behavior first, then implement until they pass. For
   pure scaffolding/config/docs tasks where a test adds no value, skip this and
   say so in `progress.log`.
4. TARGETED TESTING: there is no vitest or jest here. Node assertions come
   from tapzero; browser suites run through tapout.
   - DO NOT run the full suite (`npm test`) for every minor change. It
     rebuilds all of `dist/` and runs both the node and browser suites.
   - Node-side files run individually: `npx tsx test/0-valid.ts`, and
     likewise `1-size.ts`, `2-perf.ts`, `3-node-min.ts`.
   - Browser-side files DO NOT. `test/browser/zip.ts` and its siblings only
     export a suite function, so bundling one directly registers zero tests
     and tapout still reports a green run -- a silent false pass. To exercise
     one suite, write a temporary entry at the repo root that calls it:

     ```sh
     printf '%s\n' \
       "import * as p from './dist/browser/index.js'" \
       "import { zipSuite } from './test/browser/zip.js'" \
       "zipSuite(p as any, 'plain')" > tmp-entry.ts
     npx esbuild ./tmp-entry.ts --bundle | npx tapout --timeout 30000
     rm tmp-entry.ts
     ```

     Check the test count in the output is non-zero. This needs a current
     `dist/` (`npm run build`), and the `--timeout 30000` must stay equal to
     the one in the `test:browser` script -- see the timeout invariant in
     `CLAUDE.md`. Delete `tmp-entry.ts` before committing.
   - Run the FULL suite (`npm test`) ONLY when you believe the task is 100%
     complete, as a final gate before committing. It is the only thing that
     covers the minified bundle, where the async APIs can corrupt output
     while every sync API stays green.
   - Run only ONE npm invocation at a time: they share `dist/`, and a
     concurrent run produces a failure that looks real and is not.
5. LINT: Run `npm run lint` after any code change and fix what it reports.
6. INTEGRITY: Only consider a task done when its real tests genuinely pass.
   NEVER delete, skip, weaken, or write tautological tests to force a green
   result, and never edit `specs/prd.json` to make the loop advance without the
   work being real. If a task can't be completed honestly this session, leave
   it `passes: false` and record what blocked you in `progress.log`.
7. DOCUMENT: Update `progress.log` with what changed and any new patterns
   discovered (see format below).
8. **MARK DONE**: In `specs/prd.json`, set `passes: true` for the completed task.
   This is how the loop advances — a task is not finished until this flag is
   flipped — so include this change in the commit.
9. **COMMIT**: Once the full suite and lint pass, commit with a descriptive
   message like `FEATURE: [TaskID] - [Description]` (or `FIX:`, `TEST:`,
   `CHORE:` as appropriate).
10. ATOMICITY: Complete exactly one task per session, then stop. Do not start a
    second task even if time remains.

# PROGRESS REPORT FORMAT

APPEND to `progress.log` (never replace -- always append):

```
## [Date/Time] - [Story ID]
- What was implemented
- Files changed
- **Learnings for future iterations:**
  - Patterns discovered (e.g. "this codebase uses X for Y")
  - Gotchas encountered (e.g. "don't forget to update Z when changing W")
  - Useful context (e.g. "the evaluation panel is in component X")
---
```

The learnings section is critical — it helps future sessions avoid repeating
mistakes and understand the codebase faster.

## Consolidate Patterns

If you discover a **reusable** pattern future sessions should know, add it to
the `## Codebase Patterns` section at the TOP of `progress.log` (create it if it
doesn't exist). Keep this section tight: it is re-read in full every session, so
only consolidate **general, reusable** knowledge here — not story-specific
details, which belong in the dated entries below.

## Update AGENTS.md Files

Before committing, check whether any edited files have learnings worth
preserving in a nearby `AGENTS.md`:

1. Identify the directories you modified.
2. Check for an existing `AGENTS.md` in those directories or their parents.
3. Add valuable, reusable knowledge such as:
   - API patterns or conventions specific to that module
   - Gotchas or non-obvious requirements
   - Dependencies between files
   - Testing approaches for that area
   - Configuration or environment requirements

**Good `AGENTS.md` additions:**
- "When modifying X, also update Y to keep them in sync."
- "This module uses pattern Z for all API calls."
- "Tests require the dev server running on PORT 3000."
- "Field names must match the template exactly."

**Do NOT add:**
- Story-specific implementation details
- Temporary debugging notes
- Information already in `progress.log`

Only update `AGENTS.md` when you have genuinely reusable knowledge that would
help future work in that directory.

# STOP CONDITION

Once every story in `specs/prd.json` that is NOT flagged `"manual": true` has
`passes: true`, output the exact string `<promise>COMPLETE</promise>` as the
final thing you say, then do no further work. Manual stories are verified by a
human through `./ralph_device.sh` and never block this signal.

Do not write that string anywhere else: not in `progress.log`, not in a commit
message, and not while explaining this stop condition. The loop reads your last
message looking for it.
