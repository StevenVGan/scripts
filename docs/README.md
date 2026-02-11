# GitHub Pages (Tutorials)

This folder is the source for GitHub Pages. When enabled, it renders as HTML at `https://<username>.github.io/<repo>/`.

## Enable GitHub Pages

1. Push the repo to GitHub.
2. **Settings** → **Pages** → **Build and deployment**
3. **Source**: Deploy from a branch
4. **Branch**: `main` (or your default) → **Folder**: `/docs`
5. **Save** — GitHub will build and publish within a few minutes.

## Verify rendering

1. After deployment, visit `https://<username>.github.io/<repo>/` (e.g. `https://digan.github.io/scripts/`).
2. The index should list the CUT&RUN tutorial; the link should load the full tutorial as HTML.
3. Build status: **Settings** → **Pages** shows the deployment status; **Actions** tab shows build logs if something fails.

## Preview locally (optional)

To preview before pushing:

```bash
cd docs
bundle init
bundle add jekyll github-pages
bundle exec jekyll serve
```

Then open http://localhost:4000 (or the URL Jekyll prints).
