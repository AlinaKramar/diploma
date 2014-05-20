(TeX-add-style-hook "diploma"
 (lambda ()
    (LaTeX-add-bibliographies
     "diploma.bib")
    (TeX-run-style-hooks
     "latex2e"
     "matmex-diploma-custom10"
     "matmex-diploma-custom"
     "")))

