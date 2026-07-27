#!/usr/bin/env python3
"""Generate the vocabulary section of docs/vocabulary.html from the headers.

TempLat's public vocabulary is regular enough to be read off the source:
namespace-scope factory functions sit at exactly two spaces of indentation
(clang-format is enforced), their operand types are spelled out in a C++20
`requires` clause over a small closed set of named concepts, and almost every
header carries a `@brief` and a `Unit test: ctest -R ...` line.  This script
turns that into the type cards, the filter chips and the operation cards of
docs/vocabulary.html, so the page cannot drift out of date as the algebra
grows.

What it cannot invent is *meaning*; that comes from `@vocab-summary` lines
written into the headers.  See tools/docgen/README.md.

Usage
    gen_vocabulary.py --write      regenerate docs/vocabulary.html in place
    gen_vocabulary.py --check      exit 1 if the committed page is stale
    gen_vocabulary.py --report     audit summary, touches nothing
    gen_vocabulary.py --selftest   assert the parser still handles known traps

Exit codes: 0 ok - 1 stale - 2 hard error - 3 selftest failure.
Python >= 3.11 (tomllib), standard library only.
"""

from __future__ import annotations

import argparse
import difflib
import fnmatch
import html
import pathlib
import re
import sys
import tomllib
from collections import defaultdict

ROOT = pathlib.Path(__file__).resolve().parents[2]
CONFIG = pathlib.Path(__file__).resolve().parent / "vocabulary.toml"
PAGE = ROOT / "docs" / "vocabulary.html"

MARKERS = ("TYPES", "CONTAINERS", "TAGS", "OPERATIONS")

# A declaration head longer than this is a runaway scan, not real code.
MAX_DECL_LINES = 12

# Operators sort before names, in this order.
OPERATOR_ORDER = ["+", "-", "*", "/", "%"]

# operator symbol -> anchor-safe word
OPERATOR_SLUG = {
    "+": "plus", "-": "minus", "*": "star", "/": "slash", "%": "percent",
    "^": "caret", "&": "amp", "|": "pipe", "~": "tilde", "!": "bang",
    "<": "lt", ">": "gt", "=": "eq", "<<": "lshift", ">>": "rshift",
    "<=": "le", ">=": "ge", "==": "eqeq", "!=": "ne", "&&": "and",
    "||": "or", "++": "inc", "--": "dec", "()": "call", "[]": "index",
    "->": "arrow",
}

# Operands that mark a compile-time constant-folding overload rather than a
# distinct operation.
FOLDING_TYPES = re.compile(r"\b(ZeroType|OneType|HalfType|NumberCollection)\b")

CPP_KEYWORDS = {
    "if", "for", "while", "switch", "return", "sizeof", "alignof", "constexpr",
    "consteval", "static_assert", "decltype", "noexcept", "throw", "catch",
    "explicit", "operator", "requires", "new", "delete", "typeid",
}

BOILERPLATE = re.compile(
    r"^(a|an)\s+(class|function|struct|namespace|type)\s+(which|that)\s*",
    re.I,
)


class HardError(Exception):
    """A condition the maintainer must resolve; never papered over."""


# ───────────────────────────────────────────────────────────── lexing ──


def mask_comments(text: str) -> tuple[str, list[tuple[int, int]]]:
    """Blank out every comment, and return the doc blocks separately.

    Returns (code, docs) where `code` is the same length as `text` with all
    comment characters replaced by spaces (newlines preserved, so offsets and
    line numbers still line up), and `docs` lists the (start, end) offsets of
    every `/** ... */` block.

    Distinguishing `/**` from `/*` is what makes commented-out code invisible:
    forwardcovariantderivative.h and covariantlaplacian.h contain nothing but
    a plain block comment, so they contribute no declarations.
    """
    n = len(text)
    out = list(text)
    docs: list[tuple[int, int]] = []
    i = 0
    while i < n:
        c = text[i]
        if c in '"\'':
            # skip a string or character literal so a '//' inside one is safe
            quote, i = c, i + 1
            while i < n and text[i] != quote:
                if text[i] == "\\":
                    i += 1
                i += 1
            i += 1
            continue
        if c == "/" and i + 1 < n:
            if text[i + 1] == "/":
                j = text.find("\n", i)
                j = n if j == -1 else j
                for k in range(i, j):
                    out[k] = " "
                i = j
                continue
            if text[i + 1] == "*":
                j = text.find("*/", i + 2)
                j = n if j == -1 else j + 2
                if text.startswith("/**", i) and not text.startswith("/**/", i):
                    docs.append((i, j))
                for k in range(i, j):
                    if out[k] != "\n":
                        out[k] = " "
                i = j
                continue
        i += 1
    return "".join(out), docs


def strip_doc_stars(block: str) -> str:
    """`/** @brief x\n * y **/` -> `@brief x\ny`."""
    body = block
    body = re.sub(r"^/\*\*+", "", body)
    body = re.sub(r"\*+/$", "", body)
    lines = [re.sub(r"^\s*\*+ ?", "", ln) for ln in body.split("\n")]
    return "\n".join(lines).strip()


RE_BRIEF = re.compile(r"@brief\s+(?P<t>.*?)(?=\n\s*\n|\n\s*@|\Z)", re.S)
RE_CTEST = re.compile(r"Unit test:\s*ctest\s+-R\s+(?P<t>\S+)")
RE_VOCAB = re.compile(
    r"^@vocab(?:-(?P<key>[a-z]+(?:\+)?))?[ \t]*(?P<val>.*?)"
    r"(?=\n@vocab|\n\s*\n|\Z)",
    re.M | re.S,
)


def parse_doc(block: str) -> dict:
    body = strip_doc_stars(block)
    doc: dict = {"vocab": {}}
    m = RE_BRIEF.search(body)
    if m:
        doc["brief"] = " ".join(m.group("t").split())
    m = RE_CTEST.search(body)
    if m:
        doc["ctest"] = m.group("t")
    for m in RE_VOCAB.finditer(body):
        key = m.group("key") or "flag"
        val = " ".join(m.group("val").split())
        if key.endswith("+"):
            doc["vocab"].setdefault(key, []).append(val)
        else:
            doc["vocab"][key] = val
    return doc


# ─────────────────────────────────────────────── declaration scanning ──


def _match_angle(s: str, start: int) -> int:
    """Index just past the '>' matching the '<' at `start`, or -1."""
    depth = 0
    i = start
    while i < len(s):
        c = s[i]
        if c == "<":
            if s.startswith("<<", i):
                i += 2
                continue
            depth += 1
        elif c == ">":
            if s.startswith(">>", i) and depth >= 2:
                depth -= 2
                i += 2
                continue
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return -1


def _match_paren(s: str, start: int) -> int:
    """Index just past the ')' matching the '(' at `start`, or -1."""
    depth = 0
    for i in range(start, len(s)):
        if s[i] == "(":
            depth += 1
        elif s[i] == ")":
            depth -= 1
            if depth == 0:
                return i + 1
    return -1


def collect_declaration(lines: list[str], start: int) -> tuple[str, int, bool]:
    """Join lines from `start` until the declaration head is complete.

    Complete means: parentheses balanced, angle brackets balanced, and a '{'
    or ';' seen at depth zero.  Joining is not defensive tidiness --
    ConstructMatrix3x3's return type occupies its own line and its *name*
    starts the next line at column 2, so a line-at-a-time regex misses it.

    Returns (joined, last_line, overran).  `overran` means the head never
    closed within the cap; the caller must then advance a single line rather
    than skipping the whole window, or a bracket-counting slip in one
    declaration would silently eat the ones that follow it.
    """
    buf: list[str] = []
    paren = angle = 0
    i = start
    while i < len(lines) and i - start < MAX_DECL_LINES:
        line = lines[i]
        buf.append(line.strip())
        j = 0
        while j < len(line):
            c = line[j]
            if c == "(":
                paren += 1
            elif c == ")":
                paren -= 1
            elif line.startswith(">>", j):
                # `Outer<Inner<T>>` closes two levels; anywhere else it is a
                # right shift and means nothing to us.  Getting this wrong
                # leaves `angle` permanently positive, and the scan then runs
                # past the end of the declaration and swallows the next one.
                if angle >= 2:
                    angle -= 2
                j += 2
                continue
            elif line.startswith(("<<", "<=", ">=", "->"), j):
                j += 2
                continue
            elif c == "<":
                angle += 1
            elif c == ">":
                angle = max(0, angle - 1)
            elif c in "{;" and paren <= 0 and angle <= 0:
                return " ".join(buf), i, False
            j += 1
        i += 1
    return " ".join(buf), start, True


def strip_template_and_requires(s: str) -> tuple[str, str, str]:
    """-> (remainder, template-parameter-list, requires-clause).

    Stripping the constraint is mandatory, not cosmetic: without it
    `requires(HasEvalMethod<R> && GetNDim::get<...>() > 0) auto LatLapl(R pR)`
    parses as a function named `get`.
    """
    tparams = ""
    req = ""
    s = s.strip()
    while s.startswith("template"):
        rest = s[len("template"):].lstrip()
        if not rest.startswith("<"):
            break
        end = _match_angle(rest, 0)
        if end < 0:
            break
        tparams = rest[1:end - 1].strip()
        s = rest[end:].lstrip()
    if s.startswith("requires"):
        rest = s[len("requires"):].lstrip()
        if rest.startswith("requires"):
            # `requires requires { ... }` -- a member constraint; not ours
            return "", tparams, ""
        if rest.startswith("("):
            end = _match_paren(rest, 0)
            if end < 0:
                return "", tparams, ""
            req = rest[1:end - 1]
            s = rest[end:].lstrip()
        else:
            # `requires Concept<T> && Other<U> decl...`
            m = re.match(
                r"((?:!?\s*[\w:]+\s*(?:<[^;{]*?>)?\s*(?:&&|\|\|)?\s*)+?)"
                r"(?=\b(?:auto|constexpr|inline|static|const|[A-Z]\w*\s*<|"
                r"[A-Za-z_]\w*\s+[A-Za-z_])"
                r")",
                rest,
            )
            if m:
                req = m.group(1)
                s = rest[m.end():].lstrip()
    return s, tparams, req


RE_TRAILING_NAME = re.compile(
    r"(?:(?P<op>operator\s*(?:<<|>>|<=|>=|==|!=|&&|\|\||\+\+|--|->|\(\)|\[\]"
    r"|[-+*/%^&|~!<>=]))|(?P<id>[A-Za-z_]\w*))\s*$"
)


def parse_function(s: str) -> tuple[str, str] | None:
    """-> (name, raw parameter text), or None if `s` is not a function head."""
    depth = 0
    open_at = -1
    i = 0
    while i < len(s):
        c = s[i]
        if s.startswith(("<<", ">>", "<=", ">=", "->"), i):
            i += 2
            continue
        if c == "<":
            depth += 1
        elif c == ">":
            depth = max(0, depth - 1)
        elif c == "(" and depth == 0:
            open_at = i
            break
        i += 1
    if open_at < 0:
        return None
    head = s[:open_at].rstrip()
    m = RE_TRAILING_NAME.search(head)
    if not m:
        return None
    name = m.group("op") or m.group("id")
    name = re.sub(r"operator\s*", "operator", name)
    close = _match_paren(s, open_at)
    params = s[open_at + 1:close - 1] if close > 0 else ""
    return name, params


def split_params(params: str) -> list[str]:
    out, depth, cur = [], 0, ""
    for c in params:
        if c in "<([{":
            depth += 1
        elif c in ">)]}":
            depth -= 1
        if c == "," and depth <= 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += c
    if cur.strip():
        out.append(cur.strip())
    return out


def param_names(params: str) -> list[str]:
    names = []
    for p in split_params(params):
        p = p.replace("...", " ")
        m = re.search(r"([A-Za-z_]\w*)\s*$", p)
        names.append(m.group(1) if m else "x")
    return names


# ─────────────────────────────────────────────────────────── the scan ──


class Overload:
    __slots__ = ("name", "kind", "path", "line", "params", "tparams", "req",
                 "doc", "tags", "tag_source", "all_negated")

    def __init__(self, name, kind, path, line, params, tparams, req, doc):
        self.name = name
        self.kind = kind
        self.path = path
        self.line = line
        self.params = params
        self.tparams = tparams
        self.req = req
        self.doc = doc
        self.tags: list[str] = []
        self.tag_source = ""
        self.all_negated = False


RE_CONCEPT_ATOM = re.compile(r"(?P<neg>!\s*)?(?P<c>[A-Z]\w*)\s*<")
RE_MACRO = re.compile(r"^#define\s+(?P<name>[A-Za-z_]\w*)\s*\((?P<params>[^)]*)\)")


class Scanner:
    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.concepts: dict[str, list[str]] = cfg["concepts"]
        self.ignore: set[str] = set(cfg["concepts_ignore"]["names"])
        self.known_tags = [t["name"] for t in cfg["tags"]]
        self.exclude = set(cfg["exclude"]["names"])
        self.dir_tags = cfg["dir_tags"]
        self.section_names: dict[str, str] = {}
        for sec, body in cfg["sections"].items():
            for nm in body["names"]:
                self.section_names[nm] = sec
        self.overloads: list[Overload] = []
        self.dropped_folding: list[str] = []
        self.dropped_excluded: list[str] = []
        self.dropped_hidden: list[str] = []
        self.unknown_concepts: dict[str, set[str]] = defaultdict(set)
        self.all_function_names: set[str] = set()
        self.file_brief: dict[str, str | None] = {}
        self.overran: list[str] = []

    # -- file selection ------------------------------------------------

    def files(self) -> list[pathlib.Path]:
        out = []
        for root in self.cfg["scan"]["roots"]:
            for p in sorted((ROOT / root).rglob("*.h")):
                rel = p.relative_to(ROOT / "include").as_posix()
                if any(fnmatch.fnmatch(rel, g)
                       for g in self.cfg["scan"]["exclude_globs"]):
                    continue
                out.append(p)
        return out

    # -- main pass -----------------------------------------------------

    def scan(self) -> None:
        for path in self.files():
            rel = path.relative_to(ROOT / "include").as_posix()
            text, code, docs = self._read(path)
            self._scan_declarations(rel, code, docs)

        # Macro files come in regardless of the declaration globs -- the
        # expression macros live in util/, which is otherwise not scanned.
        for rel in self.cfg["macros"]["files"]:
            path = ROOT / "include" / rel
            if not path.exists():
                raise HardError(f"[macros].files lists a missing file: {rel}")
            text, _code, docs = self._read(path)
            self._scan_macros(rel, text, docs)

        self._drop_negative_twins()

    def _read(self, path: pathlib.Path):
        text = path.read_text(encoding="utf-8", errors="replace")
        code, blocks = mask_comments(text)
        docs = {}
        for a, b in blocks:
            docs[text.count("\n", 0, b)] = parse_doc(text[a:b])
        # A file documents one operation, so its first @brief is a usable
        # last-resort description for a free factory that carries none of
        # its own -- the docs traditionally sit on the wrapper class above.
        first = next((docs[k] for k in sorted(docs)
                      if docs[k].get("brief")), None)
        self.file_brief[path.relative_to(ROOT / "include").as_posix()] = (
            first.get("brief") if first else None
        )
        return text, code, docs

    def _attached_doc(self, doc_by_endline, line_idx: int) -> dict:
        for back in (1, 2, 3):
            d = doc_by_endline.get(line_idx - back)
            if d is not None:
                return d
        return {"vocab": {}}

    def _scan_declarations(self, rel, code, doc_by_endline) -> None:
        lines = code.split("\n")
        i = 0
        while i < len(lines):
            line = lines[i]
            if not re.match(r"^  [^\s})#*/]", line):
                i += 1
                continue
            if re.match(r"^  (namespace|using namespace|friend|return|public|"
                        r"private|protected)\b", line):
                i += 1
                continue
            # a macro body is not namespace scope, however it is indented
            if i and lines[i - 1].rstrip().endswith("\\"):
                i += 1
                continue
            joined, last, overran = collect_declaration(lines, i)
            if overran:
                self.overran.append(f"{rel}:{i + 1}")
                i += 1
                continue
            self._classify(rel, i + 1, joined,
                           self._attached_doc(doc_by_endline, i))
            i = last + 1

    def _classify(self, rel, line, joined, doc) -> None:
        s, tparams, req = strip_template_and_requires(joined)
        if not s:
            return

        m = re.match(r"(?:class|struct)\s+(\w+)\b", s)
        if m:
            self._add(rel, line, m.group(1), "type", tparams, req, doc)
            return
        m = re.match(r"using\s+(\w+)\s*=", s)
        if m:
            self._add(rel, line, m.group(1), "alias", tparams, req, doc)
            return
        if re.match(r"concept\s+\w+", s):
            return

        parsed = parse_function(s)
        if not parsed:
            return
        name, params = parsed
        if name.endswith("_impl") or name.startswith("_"):
            return
        if name in CPP_KEYWORDS:
            return  # e.g. an `if constexpr (...)` line inside a macro body
        if "= delete" in s:
            return
        self.all_function_names.add(name)
        if name in self.exclude:
            self.dropped_excluded.append(f"{name}  ({rel}:{line})")
            return
        # constant-folding overloads are implementation, not vocabulary
        if FOLDING_TYPES.search(params):
            self.dropped_folding.append(f"{name}  ({rel}:{line})")
            return
        atoms = self._concept_atoms(req, tparams)
        all_negated = bool(atoms) and all(neg for neg, _ in atoms)
        self._add(rel, line, name, "function", tparams, req, doc,
                  params=params, all_negated=all_negated)

    def _scan_macros(self, rel, text, doc_by_endline) -> None:
        allowed = set(self.cfg["macros"]["names"])
        for idx, line in enumerate(text.split("\n")):
            m = RE_MACRO.match(line)
            if not m or m.group("name") not in allowed:
                continue
            self._add(rel, idx + 1, m.group("name"), "macro", "", "",
                      self._attached_doc(doc_by_endline, idx),
                      params=m.group("params"))

    def _add(self, rel, line, name, kind, tparams, req, doc, params="",
             all_negated=False) -> None:
        if kind in ("type", "alias") and name not in self.section_names:
            return  # hundreds of expression-node classes are not vocabulary
        if "hide" in doc.get("vocab", {}):
            self.dropped_hidden.append(f"{name}  ({rel}:{line})")
            return
        ov = Overload(name, kind, rel, line, params or tparams, tparams, req, doc)
        ov.all_negated = all_negated
        self._tag(ov)
        self.overloads.append(ov)

    def _drop_negative_twins(self) -> None:
        """Remove the `requires(!Concept<T> ...)` half of a positive/negative
        overload pair.

        The idiom is a constraint and its complement in one file, where the
        complementary overload folds to ZeroType -- `LatLapl` is the canonical
        case.  The negative half is implementation, not vocabulary.

        The pairing is what makes this safe.  An all-negated clause on its own
        is a perfectly ordinary constraint: `tanh` is spelled
        `requires(!is_arithmetic_v<T> && !IsComplexType<T>)` and is the only
        overload there is, so dropping it on the negation alone would silently
        lose the operation.
        """
        positive = {(o.path, o.name) for o in self.overloads if not o.all_negated}
        keep = []
        for o in self.overloads:
            if o.all_negated and (o.path, o.name) in positive:
                self.dropped_folding.append(
                    f"{o.name}  ({o.path}:{o.line})  [complement of a "
                    f"constrained overload in the same file]")
            else:
                keep.append(o)
        self.overloads = keep

    # -- tags ----------------------------------------------------------

    def _concept_atoms(self, req: str, tparams: str = "") -> list[tuple[bool, str]]:
        """Concept applications in a requires clause, as (is_negated, name).

        The trailing '<' must actually close, otherwise `requires(M < N)` --
        a plain value comparison over two int template parameters -- reads as
        a concept named M.  Template parameters of the declaration itself are
        excluded for the same reason.
        """
        declared = set(re.findall(r"\b([A-Z]\w*)\b(?=\s*(?:=|,|\.\.\.|$))",
                                  tparams or ""))
        out = []
        for m in RE_CONCEPT_ATOM.finditer(req or ""):
            lt = req.index("<", m.start("c"))
            if _match_angle(req, lt) < 0:
                continue
            if m.group("c") in declared:
                continue
            out.append((bool(m.group("neg")), m.group("c")))
        return out

    def _tag(self, ov: Overload) -> None:
        vocab = ov.doc.get("vocab", {})
        derived: list[str] = []
        for neg, concept in self._concept_atoms(ov.req, ov.tparams):
            if concept in self.ignore:
                continue
            if concept not in self.concepts:
                self.unknown_concepts[concept].add(f"{ov.path}:{ov.line}")
                continue
            if neg:
                continue
            for t in self.concepts[concept]:
                if t not in derived:
                    derived.append(t)
        source = "requires"
        if not derived:
            derived = list(self._dir_tags(ov.path))
            source = "directory"
        if "tags" in vocab:
            derived = [t.strip() for t in vocab["tags"].split(",") if t.strip()]
            source = "@vocab-tags"
        for extra in vocab.get("tags+", []):
            for t in (x.strip() for x in extra.split(",")):
                if t and t not in derived:
                    derived.append(t)
        ov.tags = derived
        ov.tag_source = source

    def _dir_tags(self, rel: str) -> list[str]:
        best, tags = "", []
        for prefix, val in self.dir_tags.items():
            if rel.startswith(prefix) and len(prefix) > len(best):
                best, tags = prefix, val
        return tags


# ───────────────────────────────────────────────────────────── merging ──


class Entry:
    def __init__(self, name: str, file_brief: dict[str, str | None] | None = None):
        self.name = name
        self.overloads: list[Overload] = []
        self.file_brief = file_brief or {}

    # -- derived views -------------------------------------------------

    @property
    def kind(self) -> str:
        return self.overloads[0].kind

    @property
    def tags(self) -> list[str]:
        out: list[str] = []
        for ov in self.overloads:
            for t in ov.tags:
                if t not in out:
                    out.append(t)
        return out

    @property
    def primary(self) -> Overload:
        explicit = [o for o in self.overloads if "primary" in o.doc.get("vocab", {})]
        if explicit:
            return explicit[0]
        stem = self.name.lower().replace("operator", "")
        for o in self.overloads:
            if pathlib.Path(o.path).stem == self.name.lower():
                return o
        for o in self.overloads:
            if stem and pathlib.Path(o.path).stem.endswith(stem):
                return o
        return min(self.overloads, key=lambda o: (len(o.path), o.path))

    @property
    def primary_is_guessed(self) -> bool:
        if any("primary" in o.doc.get("vocab", {}) for o in self.overloads):
            return False
        return pathlib.Path(self.primary.path).stem != self.name.lower()

    def summary(self) -> tuple[str, bool]:
        """-> (text, is_authored).

        `@vocab-summary` wins.  Failing that, the operation's own `@brief`;
        failing that, the `@brief` of the file it lives in, which by the
        one-operation-per-file convention describes the same thing (the docs
        traditionally sit on the wrapper class, not the factory below it).
        """
        for o in [self.primary] + self.overloads:
            v = o.doc.get("vocab", {})
            if v.get("summary"):
                return v["summary"], True
        for o in [self.primary] + self.overloads:
            if o.doc.get("brief"):
                return _clean_brief(o.doc["brief"]), False
        for o in [self.primary] + self.overloads:
            brief = self.file_brief.get(o.path)
            if brief:
                return _clean_brief(brief), False
        return "", False

    def signatures(self) -> list[str]:
        for o in [self.primary] + self.overloads:
            v = o.doc.get("vocab", {})
            if v.get("signature"):
                return [v["signature"]]
        seen, out = set(), []
        for o in sorted(self.overloads, key=lambda o: (o.path, o.line)):
            if o.kind in ("type", "alias"):
                # `size_t _NDim = 0` is the parameter NDim, not the value 0
                tp = []
                for p in split_params(o.params):
                    p = p.split("=")[0].strip()
                    m = re.search(r"([A-Za-z_]\w*)\s*$", p)
                    if m:
                        tp.append(m.group(1).lstrip("_"))
                out.append(f"{self.name}<{', '.join(tp)}>" if tp else self.name)
                return out[:1]
            names = param_names(o.params)
            key = len(names)
            if key in seen:
                continue
            seen.add(key)
            if self.name.startswith("operator"):
                sym = self.name[len("operator"):]
                out.append(f"a {sym} b" if len(names) == 2 else f"{sym}a")
            else:
                out.append(f"{self.name}({', '.join(names)})")
        return out

    def tests(self) -> list[str]:
        out = []
        for o in self.overloads:
            t = o.doc.get("ctest")
            if t and t not in out:
                out.append(t)
        return out

    def sources(self) -> list[tuple[str, int]]:
        seen, out = set(), []
        for o in sorted(self.overloads, key=lambda o: (o.path, o.line)):
            if o.path in seen:
                continue
            seen.add(o.path)
            out.append((o.path, o.line))
        return out

    def example(self) -> str:
        for o in [self.primary] + self.overloads:
            v = o.doc.get("vocab", {})
            if v.get("example"):
                return v["example"]
        return ""


def _clean_brief(brief: str) -> str:
    """Strip the `A class which ...` header-template boilerplate off a @brief.

    Only the copula goes with it: `A class which is a classical field` should
    become `A classical field`, but `A class which implements the laplacian`
    reads better as `Implements the laplacian` than as `The laplacian`.
    """
    cleaned = BOILERPLATE.sub("", brief).strip()
    cleaned = re.sub(r"^(is|are)\s+", "", cleaned, flags=re.I)
    if not cleaned:
        return ""
    return cleaned[0].upper() + cleaned[1:]


def merge(overloads: list[Overload], split: set[str],
          file_brief: dict[str, str | None]) -> dict[str, Entry]:
    entries: dict[str, Entry] = {}
    for ov in overloads:
        key = ov.name if ov.name not in split else f"{ov.name}@{ov.path}"
        entries.setdefault(key, Entry(ov.name, file_brief)).overloads.append(ov)
    return entries


def sort_key(entry: Entry):
    if entry.name.startswith("operator"):
        sym = entry.name[len("operator"):]
        rank = OPERATOR_ORDER.index(sym) if sym in OPERATOR_ORDER else 99
        return (0, rank, sym)
    return (1, 0, entry.name.lower())


def slug(name: str) -> str:
    if name.startswith("operator"):
        sym = name[len("operator"):]
        return "op-operator-" + OPERATOR_SLUG.get(sym, "sym")
    return "op-" + re.sub(r"[^\w.-]", "-", name)


# ─────────────────────────────────────────────────────────── rendering ──


def e(s: str) -> str:
    return html.escape(s, quote=True)


def render_snip(text: str) -> str:
    """Highlight a signature the way index.html does: identifiers that are
    clearly arguments in accent, string literals in green."""
    out = e(text)
    out = re.sub(r"&quot;.*?&quot;", lambda m: f'<span class="s">{m.group(0)}</span>', out)
    return out


def render_tag_buttons(tags: list[str]) -> str:
    return "".join(
        f'<button type="button" class="tag" data-tag="{e(t)}">{e(t)}</button>'
        for t in tags
    )


def render_sources(entry: Entry, blob: str) -> str:
    bits = []
    for path, line in entry.sources():
        name = pathlib.Path(path).name
        bits.append(f'<a href="{e(blob)}include/{e(path)}#L{line}">{e(name)}</a>')
    for t in entry.tests():
        bits.append(f"<code>ctest -R {e(t)}</code>")
    return " · ".join(bits)


def render_card(entry: Entry, blob: str, extra_class: str = "") -> str:
    summary, _ = entry.summary()
    sigs = entry.signatures()
    tags = entry.tags
    cls = f"card vocab{(' ' + extra_class) if extra_class else ''}"
    parts = [
        f'<article class="{cls}" id="{slug(entry.name)}"'
        f' data-name="{e(entry.name.lower())}"'
        f' data-tags="{e(" ".join(tags))}">',
        f'  <h3><a href="#{slug(entry.name)}">{e(entry.name)}</a></h3>',
    ]
    if sigs:
        parts.append('  <pre class="snip">' + "\n".join(render_snip(s) for s in sigs) + "</pre>")
    if summary:
        parts.append(f"  <p>{e(summary)}</p>")
    if entry.example():
        parts.append(f'  <pre class="snip">{render_snip(entry.example())}</pre>')
    if tags:
        parts.append(f'  <p class="vocab-meta">{render_tag_buttons(tags)}</p>')
    src = render_sources(entry, blob)
    if src:
        parts.append(f'  <p class="vocab-src">{src}</p>')
    parts.append("</article>")
    return "\n".join(parts)


def render_chips(entries: list[Entry], tags_cfg: list[dict]) -> str:
    counts: dict[str, int] = defaultdict(int)
    for en in entries:
        for t in en.tags:
            counts[t] += 1
    out = [
        '<button type="button" class="chip" id="vocabClear" data-tag="">'
        "All <span class=\"n\">%d</span></button>" % len(entries)
    ]
    for t in tags_cfg:
        n = counts.get(t["name"], 0)
        if not n:
            continue
        out.append(
            f'<button type="button" class="chip" data-tag="{e(t["name"])}"'
            f' aria-pressed="false" title="{e(t["hint"])}">'
            f'{e(t["name"])} <span class="n">{n}</span></button>'
        )
    return "\n".join(out)


def render_section(entries: list[Entry], blob: str) -> str:
    return "\n".join(render_card(en, blob) for en in entries)


# ─────────────────────────────────────────────────────── page assembly ──


def inject(page: str, blocks: dict[str, str]) -> str:
    for key, body in blocks.items():
        pat = re.compile(
            r"(<!-- GENERATED: %s -->)(.*?)(<!-- /GENERATED: %s -->)" % (key, key),
            re.S,
        )
        if not pat.search(page):
            raise HardError(
                f"docs/vocabulary.html has no <!-- GENERATED: {key} --> "
                f"... <!-- /GENERATED: {key} --> marker pair"
            )
        page = pat.sub(lambda m: m.group(1) + "\n" + body + "\n" + m.group(3),
                       page, count=1)
    return page


def build(cfg: dict) -> tuple[dict[str, str], Scanner, dict[str, Entry]]:
    scanner = Scanner(cfg)
    scanner.scan()

    split = set(cfg["split"]["names"])
    entries = merge(scanner.overloads, split, scanner.file_brief)

    # -- hard errors ---------------------------------------------------
    missing = []
    for sec, body in cfg["sections"].items():
        for nm in body["names"]:
            if nm not in entries:
                missing.append(f"{nm}  (listed in [sections.{sec}])")
    if missing:
        raise HardError(
            "these names are configured but no longer found in the headers "
            "-- they were renamed, deleted, or the scan globs changed:\n  "
            + "\n  ".join(missing)
        )

    known = {t["name"] for t in cfg["tags"]}
    bad = []
    for en in entries.values():
        for t in en.tags:
            if t not in known:
                bad.append(f"{t!r} on {en.name} ({en.primary.path})")
    if bad:
        raise HardError(
            "unknown tag(s) -- add them to [[tags]] or fix the "
            "@vocab-tags line:\n  " + "\n  ".join(sorted(set(bad)))
        )

    if scanner.unknown_concepts:
        lines = [
            f"{c}  (e.g. {sorted(w)[0]})"
            for c, w in sorted(scanner.unknown_concepts.items())
        ]
        raise HardError(
            "new concept(s) appear in a requires clause but are in neither "
            "[concepts] nor [concepts_ignore].\nDecide whether each one means "
            "a new filter tag, then add it to one of the two:\n  "
            + "\n  ".join(lines)
        )

    # -- section partition ---------------------------------------------
    claimed = {sec: [entries[n] for n in cfg["sections"][sec]["names"]]
               for sec in cfg["sections"]}
    in_a_card_section = {
        n for sec in ("types", "collections", "factories")
        for n in cfg["sections"][sec]["names"]
    }
    ops = sorted(
        (en for k, en in entries.items() if en.name not in in_a_card_section),
        key=sort_key,
    )

    blob = cfg["github"]["blob_base"]
    blocks = {
        "TYPES": render_section(claimed["types"], blob),
        "CONTAINERS": (
            render_section(claimed["collections"], blob)
            + '\n</div>\n\n<h3 class="vocab-subhead">Constructors</h3>\n'
            + '<div class="feat-grid">\n'
            + render_section(claimed["factories"], blob)
        ),
        "TAGS": render_chips(ops, cfg["tags"]),
        "OPERATIONS": render_section(ops, blob),
    }
    return blocks, scanner, entries


# ──────────────────────────────────────────────────────────── reporting ──


def report(scanner: Scanner, entries: dict[str, Entry], cfg: dict) -> None:
    ops = [en for en in entries.values()]
    print(f"headers scanned          {len(scanner.files())}")
    print(f"overloads kept           {len(scanner.overloads)}")
    print(f"entries after merge      {len(entries)}")
    print(f"dropped: folding         {len(scanner.dropped_folding)}")
    print(f"dropped: [exclude]       {len(scanner.dropped_excluded)}")
    print(f"dropped: @vocab-hide     {len(scanner.dropped_hidden)}")
    print(f"distinct function names  {len(scanner.all_function_names)}")
    print(f"declaration heads unclosed {len(scanner.overran)}")

    no_summary = sorted(en.name for en in ops if not en.summary()[0])
    derived_only = sorted(en.name for en in ops
                          if en.summary()[0] and not en.summary()[1])
    dir_tagged = sorted({ov.name for ov in scanner.overloads
                         if ov.tag_source == "directory"})
    guessed = sorted(en.name for en in ops
                     if len(en.overloads) > 1 and en.primary_is_guessed)

    def block(title, items):
        print(f"\n{title}  ({len(items)})")
        for it in items:
            print("   ", it)

    block("no summary at all -- card will render bare", no_summary)
    block("summary derived from @brief -- consider @vocab-summary", derived_only)
    block("tags fell back to the directory (no requires clause)", dir_tagged)
    block("primary overload guessed -- consider @vocab-primary", guessed)
    block("dropped as constant folding", sorted(scanner.dropped_folding))
    block("dropped by [exclude].names", sorted(scanner.dropped_excluded))


# ───────────────────────────────────────────────────────────── selftest ──


def selftest(entries: dict[str, Entry], scanner: Scanner) -> int:
    fails = []

    def ok(cond, msg):
        if not cond:
            fails.append(msg)

    # Every declaration head must close. One that does not means the bracket
    # counting slipped, and the scan then eats the declarations that follow --
    # which is how asComplexField, dag, getAverager and projectRadiallyFourier
    # silently went missing when `>>` was treated as a right shift.
    ok(not scanner.overran,
       f"{len(scanner.overran)} declaration head(s) never closed, e.g. "
       f"{scanner.overran[:3]}")

    # multi-line return type: the name starts its own line at column 2
    ok("ConstructMatrix3x3" in entries,
       "ConstructMatrix3x3 missing -- multi-line declaration join is broken")

    # the four that the `>>` bug used to swallow
    for swallowed in ("asComplexField", "dag", "projectRadiallyFourier"):
        ok(swallowed in entries,
           f"{swallowed} missing -- a preceding declaration is being over-read")

    # macro pass + allowlist
    for m in ("MakeVector", "Total"):
        ok(m in entries, f"{m} missing -- macro scan is broken")

    # entirely commented-out files must contribute nothing
    for dead in ("ForwardCovariantDerivative", "CovariantLaplacian",
                 "ComplexLaplacian"):
        ok(dead not in entries,
           f"{dead} present -- commented-out code is being parsed")

    # unconstrained factory + cross-directory merge
    if "trace" in entries:
        t = set(entries["trace"].tags)
        ok("SU2Field" in t,
           f"trace tags {sorted(t)} -- expected the su2algebra directory fallback")
    else:
        fails.append("trace missing")

    # ConditionalBinaryGetter is defined exclusively (it negates every other
    # Has*Get), so the generic operator+ must tag as plain Field and nothing
    # else -- and add.h's four ZeroType/HalfType folding overloads must be
    # gone, leaving exactly one kept overload in that file.
    add = [o for o in scanner.overloads
           if o.name == "operator+"
           and o.path == "TempLat/lattice/algebra/operators/add.h"]
    ok(len(add) == 1,
       f"operators/add.h kept {len(add)} operator+ overloads -- expected 1 "
       f"(the folding overloads should be dropped)")
    if add:
        ok(add[0].tags == ["Field"],
           f"generic operator+ tags {add[0].tags} -- expected exactly ['Field']")

    # ...while the merged card unions every algebra that defines the symbol.
    if "operator+" in entries:
        t = set(entries["operator+"].tags)
        ok({"Field", "ComplexField", "SU2Field"} <= t,
           f"merged operator+ tags {sorted(t)} -- expected at least "
           f"Field, ComplexField and SU2Field")
    else:
        fails.append("operator+ missing")

    # merge policy: shift spans five algebras
    if "shift" in entries:
        n = len(entries["shift"].tags)
        ok(n >= 5, f"shift has {n} tags -- expected >= 5 (merge across algebras)")
    else:
        fails.append("shift missing")

    for f in fails:
        print("FAIL:", f, file=sys.stderr)
    if not fails:
        print("selftest: all assertions passed")
    return 3 if fails else 0


# ───────────────────────────────────────────────────────────────── main ──


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", action="store_true", help="regenerate the page")
    g.add_argument("--check", action="store_true", help="fail if page is stale")
    g.add_argument("--report", action="store_true", help="audit summary only")
    g.add_argument("--selftest", action="store_true", help="parser assertions")
    args = ap.parse_args()

    cfg = tomllib.loads(CONFIG.read_text(encoding="utf-8"))

    try:
        blocks, scanner, entries = build(cfg)
    except HardError as exc:
        print(f"\ngen_vocabulary: {exc}\n", file=sys.stderr)
        return 2

    if args.report:
        report(scanner, entries, cfg)
        return 0
    if args.selftest:
        return selftest(entries, scanner)

    if not PAGE.exists():
        print(f"gen_vocabulary: {PAGE} does not exist", file=sys.stderr)
        return 2
    current = PAGE.read_text(encoding="utf-8")
    try:
        updated = inject(current, blocks)
    except HardError as exc:
        print(f"\ngen_vocabulary: {exc}\n", file=sys.stderr)
        return 2

    if args.write:
        if updated != current:
            PAGE.write_text(updated, encoding="utf-8")
            print(f"gen_vocabulary: wrote {PAGE.relative_to(ROOT)}")
        else:
            print("gen_vocabulary: already up to date")
        return 0

    if updated == current:
        print("gen_vocabulary: docs/vocabulary.html is up to date")
        return 0
    diff = difflib.unified_diff(
        current.splitlines(True), updated.splitlines(True),
        "docs/vocabulary.html (committed)", "docs/vocabulary.html (regenerated)",
    )
    sys.stderr.writelines(diff)
    print("\ngen_vocabulary: docs/vocabulary.html is STALE -- run --write",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
