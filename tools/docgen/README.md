# docgen — the vocabulary page generator

`docs/vocabulary.html` lists every field type, container and expression-algebra
operation TempLat provides, filterable by the types that take part in each one.
Most of it is read out of the headers by `gen_vocabulary.py`, so it does not rot
as the library grows.

```sh
python3 tools/docgen/gen_vocabulary.py --write      # regenerate the page
python3 tools/docgen/gen_vocabulary.py --check      # exit 1 if it is stale
python3 tools/docgen/gen_vocabulary.py --report     # audit summary
python3 tools/docgen/gen_vocabulary.py --selftest   # parser assertions
```

Python ≥ 3.11, standard library only. Run `--write` and commit the result
whenever you add or change an operation — same habit as `./format.sh`. CI
(`.github/workflows/docs-vocabulary.yml`) fails the PR if you forget.

## What is automatic

Adding a new operation costs **no annotation at all** to appear on the page
correctly tagged. The generator derives:

| | from |
|---|---|
| the name and signature | the factory declaration |
| **the filter tags** | the concept names in its `requires` clause |
| the description | `@brief`, falling back to the file's `@brief` |
| the source link | the file and line |
| the test to run | the `Unit test: ctest -R …` line |

Tags are the interesting one. Operand types are spelled out as C++20 concepts —
`HasSU2Get`, `HasComplexFieldGet`, `ConditionalListBinaryGetter` — and
`vocabulary.toml` maps each to a tag. Negated atoms (`!HasSU2Get<T>`, present
only to break overload ambiguity) are ignored. A factory with no constraint
falls back to a tag derived from its directory.

Three things make the page **fail loudly** rather than quietly degrade:

1. a name listed in `[sections.*]` that no longer exists — catches renames;
2. a `@vocab-tags` value that is not a declared tag — catches typos;
3. **a concept in a `requires` clause that is in neither `[concepts]` nor
   `[concepts_ignore]`** — a new concept forces an explicit decision about
   whether it means a new filter tag.

## What you write by hand

Meaning. Put it in the doc block on the **factory function**, not the wrapper
class above it, so that whoever changes the `requires` clause sees the
description in the same diff.

```cpp
/** @brief A class which implements the laplacian.
 *
 * @vocab-summary Lattice Laplacian: the standard $2d+1$-point stencil, using the operand's own $dx$.
 * @vocab-signature LatLapl(expr)
 * @vocab-example rho = 0.5 * pow<2>(pi) - 0.5 * phi * LatLapl(phi);
 * @vocab-tags+ Fourier
 * @vocab-primary
 *
 * Unit test: ctest -R test-laplacianlocal
 **/
```

| tag | effect |
|---|---|
| `@vocab-summary` | the card's prose. `$…$` renders as maths via KaTeX. |
| `@vocab-signature` | replaces the synthesised signature. May span lines. |
| `@vocab-tags` | **replaces** the derived tag set, comma-separated. |
| `@vocab-tags+` | **adds** to it. |
| `@vocab-example` | a second code line under the summary. |
| `@vocab-primary` | this overload's doc wins when several files define the name. |
| `@vocab-hide` | omit the entry entirely. |
| `@vocab-name` | override the displayed name. |

A `@vocab-*` tag on a **member** function makes it eligible despite its
indentation — the escape hatch for things like `inFourierSpace()`.

## Structure

`vocabulary.toml` holds the *page* decisions: which types belong in which
section and in what order, the concept→tag dictionary, exclusions, and the
macro allowlist. In-source tags describe the thing; the config describes the
page.

Classes are rendered **only** if named in `[sections.*]` — the scan sees several
hundred expression-node classes that are implementation, not vocabulary.

`docs/vocabulary.html` is hand-written apart from four regions:

```html
<!-- GENERATED: TYPES -->  <!-- GENERATED: CONTAINERS -->
<!-- GENERATED: TAGS -->   <!-- GENERATED: OPERATIONS -->
```

Only their interiors are replaced, so prose and layout stay freely editable and
a diff shows exactly which cards changed. Output is static HTML rather than a
JSON payload, so the page works over `file://` and without JavaScript;
filtering is progressive enhancement layered on the DOM that is already there.

## Known limits

- Member functions are invisible unless annotated (the two-space indent rule is
  what makes the scanner reliable).
- Whether an operation is meaningful in Fourier space is a physics fact absent
  from the type system: `@vocab-tags+ Fourier`, by hand.
- Section 3 is sorted alphabetically because deterministic output is what makes
  `--check` usable. "The five things you need first" belongs in the hand-written
  `<p class="lede">` above the list.
