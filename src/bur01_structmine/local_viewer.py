"""Embed the pinned 3Dmol bundle without a CDN or asynchronous download."""
from functools import lru_cache

from .config import project_root


@lru_cache(maxsize=1)
def _bundle() -> str:
    return (project_root() / "ui" / "assets" / "3Dmol-2.5.5-min.js").read_text(encoding="utf-8")


def offline_viewer_html(viewer) -> str:
    # Resolve the py3Dmol loader before its guarded asynchronous CDN branch.
    bundle = _bundle().replace("</script", "<\\/script")
    return ('<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<style>body{margin:0}canvas{max-width:100%}</style><script>' + bundle +
            '</script><script>var $3Dmolpromise=Promise.resolve();</script>' + viewer._make_html())
