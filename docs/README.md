# GitHub Pages (Tutorials)

This folder is the source for GitHub Pages. Uses the **[Minima](https://github.com/jekyll/minima)** theme — Jekyll's default theme (simple, minimal).

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
   layout: page
   title: Your Tutorial Title
   permalink: /tutorials/YourTutorial/
   ---
   ```
2. Add to `_config.yml` under `minima.nav_pages` if you want it in the header:
   ```yaml
   minima:
     nav_pages:
       - tutorials/YourTutorial.md
   ```

## Preview locally

```bash
cd docs
bundle init
bundle add jekyll jekyll-feed jekyll-seo-tag
bundle exec jekyll serve
```

Then open http://localhost:4000.
