(TeX-add-style-hook "presentation"
 (lambda ()
    (TeX-add-symbols
     "genmap")
    (TeX-run-style-hooks
     "booktabs"
     "listings"
     "graphicx"
     "polyglossia"
     "hyperref"
     "xecyr"
     "xltxtra"
     "caption"
     "xunicode"
     "fontspec"
     "latex2e"
     "beamer10"
     "beamer"
     "")))

