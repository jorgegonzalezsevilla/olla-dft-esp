# Galería de plantillas visuales

La misma figura de bandas y DOS del silicio en cada estilo (`templates`,
`plot -t`).

    olla-dft templates list              # ver las disponibles
    olla-dft templates show dark         # qué define cada una
    olla-dft plot calc/ -t latex-true    # usarla

Para una plantilla propia:

    olla-dft templates export dark       # escribe un JSON editable
    olla-dft plot calc/ -t dark-copia    # la usa por nombre

### Archivos

| Archivo | Qué es |
|---|---|
| `galeria_plantillas.png` | comparación de todas las plantillas |
| `latex.pdf` | Computer Modern sin necesidad de LaTeX instalado |
| `latex-true.pdf` | renderizado con LaTeX real (`--usetex`) |
| `dark.pdf` | tema oscuro para diapositivas |
