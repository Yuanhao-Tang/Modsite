Please read `LEGACY_SOURCES.md` first.

This project is a formal R package rebuild.
The files listed in `LEGACY_SOURCES.md` are reference materials for statistical logic,
domain behavior, legacy package structure, and migration context only.
They are not the target architecture for the new package.

Requirements:
- Use English only.
- Use formal R package conventions.
- Use a consistent snake_case API unless there is a strong object-system reason not to.
- Do not copy legacy files verbatim unless explicitly requested.
- Preserve validated statistical intent, but feel free to redesign naming, layout, and module boundaries.

Before editing, first produce:
1. a proposed file organization plan
2. an old-to-new function rename map
3. a phased migration plan
