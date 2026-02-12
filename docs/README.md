# GitHub Pages (Tutorials)

This folder is the source for GitHub Pages. Uses the **[Minima](https://github.com/jekyll/minima)** theme — Jekyll's default theme (simple, minimal).

## Enable GitHub Pages

1. Push the repo to GitHub.
2. **Settings** → **Pages** → **Build and deployment**
3. **Source**: Deploy from a branch
4. **Branch**: `main` (or your default) → **Folder**: `/docs`
5. **Save** — GitHub will build and publish within a few minutes.

## Adding new tutorials

Each tutorial lives in its own subdirectory under `tutorials/` so files stay organized:

1. Create `tutorials/your_tutorial/index.md` with front matter:
   ```yaml
   ---
   layout: page
   title: Your Tutorial Title
   permalink: /tutorials/your_tutorial/
   ---
   ```
2. Place tutorial-specific assets (images, data files) in the same folder.
3. Add to `_config.yml` under `minima.nav_pages` if you want it in the header:
   ```yaml
   minima:
     nav_pages:
       - tutorials/your_tutorial/index.md
   ```
4. Add the tutorial to the index page:
   ```markdown
   ## Available tutorials

   - **[Your Tutorial](tutorials/your_tutorial/)** — Your tutorial description.
   ```
