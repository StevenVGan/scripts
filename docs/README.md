# GitHub Pages (Tutorials)

This folder is the source for GitHub Pages. Uses the **[Just the Docs](https://just-the-docs.github.io/just-the-docs/)** theme — a Jekyll documentation theme with built-in sidebar, TOC, collapsible nav, and code highlighting.

## Enable GitHub Pages

1. Push the repo to GitHub.
2. **Settings** → **Pages** → **Build and deployment**
3. **Source**: Deploy from a branch
4. **Branch**: `main` (or your default) → **Folder**: `/docs`
5. **Save** — GitHub will build and publish within a few minutes.

## Adding new tutorials

1. Create `tutorials/YourTutorial.md` with front matter:
   ```yaml
   ---
   layout: default
   title: Your Tutorial Title
   ---
   ```
2. Add to `_config.yml` under `nav`:
   ```yaml
   - title: Your Tutorial
     url: /tutorials/YourTutorial/
   ```

## Preview locally

```bash
cd docs
bundle init
bundle add jekyll jekyll-remote-theme just-the-docs
bundle exec jekyll serve
```

Then open http://localhost:4000.
